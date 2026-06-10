using Oceananigans
using JLD2
using Printf
using Statistics

using Oceananigans.Units: hours, hour, days
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using Oceananigans.TurbulenceClosures

################################################
#### Surface wind boundary condition ####
################################################

@inline function surface_u10_profile(y, p)
    lat = y / 111e3

    northern_trade = exp(-((lat - p.jet_lat) / p.jet_width)^2)
    southern_trade = exp(-((lat + p.jet_lat) / p.jet_width)^2)

    trade_shape = max(northern_trade, southern_trade)

    return p.u₁₀_equator + (p.u₁₀_trades - p.u₁₀_equator) * trade_shape
end

@inline function spinup_ramp(t, p)
    return min(1.0, t / p.spinup_time)
end

# Resting-ocean ("ocean at rest") bulk stress — kept as a reference profile for
# plotting; the model now uses the relative-wind stress below.
@inline function τx_profile(x, y, t, p)
    u10 = surface_u10_profile(y, p)
    τ = p.ρₐ / p.ρₒ * p.cᴰ * u10 * abs(u10)
    return τ * spinup_ramp(t, p)
end

################################################
#### Relative-wind (current-feedback) stress ####
################################################
# Full bulk stress using the wind RELATIVE to the ocean surface current:
#   τ = (ρₐ/ρₒ) C_D (U_air − U_ocean) |U_air − U_ocean|.
# Sign bookkeeping: work in an eastward-positive frame, where the easterly
# trades are NEGATIVE, so the atmospheric zonal wind is −(wind speed).
# Oceananigans' top momentum flux has the OPPOSITE sign to the physical stress
# (a positive top u-flux drives WESTWARD flow), hence the leading minus on the
# return.  At rest (u = v = 0) this reduces exactly to τx_profile above.
#
# The BCs are field-dependent on (u, v): Oceananigans interpolates both to each
# component's boundary location automatically (same mechanism as the T and S
# fluxes, which depend on T and S).
@inline atmospheric_u10(y, p) = -surface_u10_profile(y, p)   # eastward-positive (trades < 0)

@inline function τx_relative(x, y, t, u, v, p)
    Δu    = atmospheric_u10(y, p) - u      # relative zonal wind
    Δv    = -v                             # relative meridional wind (no atmospheric v)
    speed = sqrt(Δu^2 + Δv^2)              # |U_air − U_ocean|
    τphys = p.ρₐ / p.ρₒ * p.cᴰ * Δu * speed
    return -τphys * spinup_ramp(t, p)
end

@inline function τy_relative(x, y, t, u, v, p)
    Δu    = atmospheric_u10(y, p) - u
    Δv    = -v
    speed = sqrt(Δu^2 + Δv^2)
    τphys = p.ρₐ / p.ρₒ * p.cᴰ * Δv * speed
    return -τphys * spinup_ramp(t, p)
end

# NOTE: the prescribed north/south wall temperature fluxes were removed — the
# N/S walls are now INSULATING (zero flux), consistent with the salinity walls
# and the bottom.  Their job (net-zero vertical heat redistribution at the
# subtropical edges) is superseded by the poleward heat sink (§ ns_band / sink_T),
# which provides the actual net heat export and sets the cool subtropical column.

################################################
#### Surface heat flux (SST restoring) ####
################################################
# Oceananigans sign convention: a POSITIVE top T-flux removes heat (cooling).
# So a stable relaxation toward T_target must be  F = +λ_T (T - T_target):
# when T > T_target the flux is positive → cooling → pulls T back down
# (NEGATIVE feedback).  Equilibrium SST → T_target.
#
# ⚠ Earlier this had the sign flipped (F = Q₀ - λ_T(T-T_restore)), which is an
# ANTI-restoring positive feedback — it made the equilibrium a repeller and the
# surface temperature ran away (hot → 156 °C crash, or cold → homogenised).
#
# The restoring target is ZONALLY STRUCTURED: a warm pool in the west and a
# cold-tongue target in the east, with a front in the central/eastern basin.
# A uniform target pins the whole surface to one value and erases the SST
# gradient; this anchors the warm-west / cold-east structure (and holds it
# against the slow warming of the eastern thermocline / upwelling source).
@inline function T_target_profile(x, p)
    xnorm = x / p.Lx
    warm_fraction = 0.5 * (1 - tanh((xnorm - p.x_front) / p.front_width))
    return p.T_cold + (p.T_warm - p.T_cold) * warm_fraction
end

# Evaporation rate E(x) [m s⁻¹], enhanced in the warm pool (west).  ONE field
# drives both the latent cooling here and the salinification in Jˢ.
@inline function evaporation_rate_profile(x, p)
    xnorm = x / p.Lx
    warm  = 0.5 * (1 - tanh((xnorm - p.x_front) / p.front_width))  # 1 west → 0 east
    return p.E_base + p.E_warm * warm
end

# Surface heat flux = SST restoring + evaporative (latent) cooling.
# Oceananigans convention: +ve flux = heat OUT.  Latent cooling is +ve (heat
# leaves the ocean as the warm pool evaporates), so it cools the surface.
@inline function surface_T_flux_bc(x, y, t, T, p)
    restoring      = p.λ_T * (T - T_target_profile(x, p))
    latent_cooling = p.latent_coeff * evaporation_rate_profile(x, p)
    solar_heating  = -p.Q_solar          # -ve flux = heat INTO the ocean
    return restoring + latent_cooling + solar_heating
end

################################################
#### Boundary sponge — damp u and v only ####
################################################

@inline function sponge_mask(x, y, p)
    west  = max(0, (p.sponge_width_x - x) / p.sponge_width_x)
    east  = max(0, (x - (p.Lx - p.sponge_width_x)) / p.sponge_width_x)
    south = max(0, (p.sponge_width_y - (y + p.Ly/2)) / p.sponge_width_y)
    north = max(0, (y - (p.Ly/2 - p.sponge_width_y)) / p.sponge_width_y)
    east_west_mask   = max(west^2, east^2)
    north_south_mask = max(south^2, north^2)
    return p.east_west_strength * east_west_mask +
           p.north_south_strength * north_south_mask
end

# Depth-dependent Rayleigh drag: ~0 above z_damp_top, ramps to 1 below
# z_damp_bottom.  Relaxes deep flow to rest so zonal momentum does not
# accumulate below the EUC and the EUC stays confined to the upper ocean.
@inline function deep_mask(z, p)
    return clamp((p.z_damp_top - z) / (p.z_damp_top - p.z_damp_bottom), 0.0, 1.0)
end

@inline sponge_u(x, y, z, t, u, p) =
    -(p.λ₀ * sponge_mask(x, y, p) + p.λ_deep * deep_mask(z, p)) * u
@inline sponge_v(x, y, z, t, v, p) =
    -(p.λ₀ * sponge_mask(x, y, p) + p.λ_deep * deep_mask(z, p)) * v

################################################
#### Poleward heat sink (N/S temperature relaxation) ####
################################################
# Relaxes T toward a cool "extratropical" reference profile within the N/S
# sponge bands — a poleward heat EXPORT.  Warm water carried to the subtropical
# edges is cooled there, removing the net tropical heat uptake (~+1.15 PW) that
# the cold-tongue restoring injects.  Without this sink that heat just
# accumulates and the basin warms indefinitely instead of equilibrating.

# N/S-only mask (1 near the north & south edges, 0 in the interior).
@inline function ns_band_mask(y, p)
    south = max(0, (p.sponge_width_y - (y + p.Ly / 2)) / p.sponge_width_y)
    north = max(0, (y - (p.Ly / 2 - p.sponge_width_y)) / p.sponge_width_y)
    return max(south^2, north^2)
end

# Cool extratropical reference column: surface T_surf_ref, deep T_deep_ref.
@inline function T_extratropical_ref(z, p)
    return p.T_deep_ref + (p.T_surf_ref - p.T_deep_ref) *
           0.5 * (1 + tanh((z - p.z_ref) / p.d_ref))
end

@inline sink_T(x, y, z, t, T, p) =
    -p.λ_T_sink * ns_band_mask(y, p) * (T - T_extratropical_ref(z, p))

################################################
#### Salinity boundary condition ####
################################################

@inline function Jˢ(x, y, t, S, p)
    evaporation_flux = -evaporation_rate_profile(x, p) * S   # -ve → salinity increases
    restoring_flux   = -p.restore_rate * (S - p.S_surface_target)
    return evaporation_flux + restoring_flux
end

################################################
#### Build model ####
################################################

function build_model()

    Nx = 140
    Ny = 40
    Nz = 60

    Lx = 1.5e7
    Ly = 2.22e6
    Lz = 2000

    refinement = 3.5
    stretching = 18

    h(k)  = (k - 1) / Nz
    ζ₀(k) = 1 + (h(k) - 1) / refinement
    Σ(k)  = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))
    z_interfaces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

    grid = RectilinearGrid(CPU(),
                           size     = (Nx, Ny, Nz),
                           halo     = (4, 4, 4),   # ≥4 required for WENO(order=7)
                           x        = (0, Lx),
                           y        = (-Ly/2, Ly/2),
                           z        = z_interfaces,
                           topology = (Bounded, Bounded, Bounded))

    ρₒ = 1026
    equation_of_state = TEOS10EquationOfState(reference_density = ρₒ)
    buoyancy = SeawaterBuoyancy(; equation_of_state)

    ################################################
    #### Temperature boundary conditions ####
    ################################################

    cᴾ  = 3991

    # SST-restoring piston velocity.  λ_T [m s⁻¹] = (λ in W m⁻² K⁻¹)/(ρₒ cᴾ).
    # λ_T = 7e-6 ≈ 29 W m⁻² K⁻¹ (a realistic air-sea damping) → restoring timescale
    # over a 50 m layer ≈ 83 days.  WEAKENED from 2e-5 (≈80 W m⁻² K⁻¹, a stiff
    # thermostat) so the ocean dynamics + evaporative cooling actually shape SST
    # instead of being overruled by the restoring.  Target SST is zonally
    # structured: warm pool (T_warm) west, cold tongue (T_cold) east, front at x_front.
    λ_T         = 7e-6
    T_warm      = 29.5     # western warm-pool target
    T_cold      = 24.0     # eastern cold-tongue target
    x_front     = 0.6      # front location as fraction of basin width (W→E)
    front_width = 0.18     # front sharpness (fraction of Lx)

    # Warm-pool evaporation → latent cooling (a physical heat sink) + salinification.
    # E(x) = E_base + E_warm·warm_fraction(x).  E_warm = 3e-8 m/s (~2.6 mm/day) gives
    # ρ_fw Lᵥ E /(ρₒ cᴾ) ≈ 100 W m⁻² of cooling in the warm pool (~25 W m⁻² in the east).
    Lᵥ           = 2.5e6      # J kg⁻¹, latent heat of vaporisation
    ρ_fw         = 1000.0     # kg m⁻³, freshwater density
    E_base       = 1e-8       # background evaporation (m s⁻¹)
    E_warm       = 3e-8       # warm-pool evaporation enhancement (m s⁻¹)
    latent_coeff = ρ_fw * Lᵥ / (ρₒ * cᴾ)   # E (m/s) → kinematic heat flux (°C m s⁻¹)

    # Constant solar heat input (heat IN), set to roughly balance the area-mean
    # evaporative cooling (~70 W m⁻²) so the net surface heat budget ≈ 0 and SST
    # stays near target.  The restoring + N/S fluxes then handle the spatial
    # redistribution.  Tune Q_solar_Wm2 against the surface_heat_budget diagnostic.
    Q_solar_Wm2 = 70.0
    Q_solar     = Q_solar_Wm2 / (ρₒ * cᴾ)  # kinematic (°C m s⁻¹)

    surface_T_flux_parameters = (; λ_T, T_warm, T_cold, x_front, front_width, Lx,
                                   E_base, E_warm, latent_coeff, Q_solar)

    surface_T_bc = FluxBoundaryCondition(surface_T_flux_bc,
                                         field_dependencies = :T,
                                         parameters         = surface_T_flux_parameters)

    T_bcs = FieldBoundaryConditions(
        top    = surface_T_bc,
        bottom = FluxBoundaryCondition(0)
        # north & south default to zero flux (insulating walls)
    )

    ################################################
    #### Salinity boundary conditions ####
    ################################################

    restore_rate     = 5e-7
    S_surface_target = 34.6

    # Evaporation (E_base, E_warm, front shape) is SHARED with the heat flux, so the
    # same warm-pool evaporation that cools the surface also makes it saltier.
    salinity_flux_parameters = (; restore_rate, S_surface_target,
                                  E_base, E_warm, x_front, front_width, Lx)

    evaporation_bc = FluxBoundaryCondition(Jˢ,
                                           field_dependencies = :S,
                                           parameters         = salinity_flux_parameters)

    S_bcs = FieldBoundaryConditions(top = evaporation_bc)

    ################################################
    #### Wind boundary condition ####
    ################################################

    u₁₀_trades  = 8
    jet_lat     = 5.0
    jet_width   = 2.5
    cᴰ          = 2e-3
    ρₐ          = 1.2
    spinup_time = 10days

    τx_parameters = (;
        u₁₀_trades,
        jet_lat,
        jet_width,
        cᴰ,
        ρₐ,
        ρₒ,
        Ly,
        spinup_time
    )

    # Relative-wind (current-feedback) stress: both components depend on the surface
    # (u, v).  τy_relative opposes the ocean's own v (eddy damping) — v now carries a
    # top stress where the resting-ocean formulation gave it none.  At rest this
    # reduces exactly to τx_profile, so it's a clean add-on over the sink-fixed base.
    u_bcs = FieldBoundaryConditions(
        top    = FluxBoundaryCondition(τx_relative, field_dependencies = (:u, :v),
                                       parameters = τx_parameters),
        bottom = FluxBoundaryCondition(0)
    )

    v_bcs = FieldBoundaryConditions(
        top    = FluxBoundaryCondition(τy_relative, field_dependencies = (:u, :v),
                                       parameters = τx_parameters),
        bottom = FluxBoundaryCondition(0)
    )

    ################################################
    #### Sponge forcing ####
    ################################################

    sponge_width_x       = 5.0e5
    sponge_width_y       = 5.0e5
    sponge_time          = 5days
    λ₀                   = 1 / sponge_time
    east_west_strength   = 1.0
    north_south_strength = 1.0

    # Deep Rayleigh drag: zero in the upper ~300 m (where the EUC lives),
    # full strength below ~550 m.  λ_deep sets the deep spin-down timescale.
    z_damp_top    = -300.0
    z_damp_bottom = -550.0
    λ_deep        = 1 / 20days

    sponge_params = (;
        sponge_width_x,
        sponge_width_y,
        Lx,
        Ly,
        λ₀,
        east_west_strength,
        north_south_strength,
        z_damp_top,
        z_damp_bottom,
        λ_deep
    )

    u_forcing = Forcing(sponge_u; field_dependencies = :u, parameters = sponge_params)
    v_forcing = Forcing(sponge_v; field_dependencies = :v, parameters = sponge_params)

    # Poleward heat sink: relax T toward a cool extratropical profile in the N/S
    # bands so the net tropical heat uptake is exported and the basin equilibrates.
    # λ_T_sink sets the sink timescale; T_surf_ref the cool subtropical surface.
    # Tune against the surface_heat_budget diagnostic (dHC/dt should → 0).
    λ_T_sink   = 1 / 150days   # weakened from 1/20days — that removed ~9 PW (≈8× the
                               # ~1.15 PW target) and over-cooled the whole basin by ~5°C.
    T_surf_ref = 22.0      # cool extratropical surface (°C)
    T_deep_ref = 3.5       # deep reference (matches initial deep T → no deep nudging)
    z_ref      = -150.0    # reference thermocline depth (m)
    d_ref      = 120.0     # reference thermocline thickness (m)

    sink_params = (;
        sponge_width_y,
        Ly,
        λ_T_sink,
        T_surf_ref,
        T_deep_ref,
        z_ref,
        d_ref
    )

    T_forcing = Forcing(sink_T; field_dependencies = :T, parameters = sink_params)

    ################################################
    #### Model ####
    ################################################

    Ω = 7.2921e-5
    R = 6.371e6
    β = 2Ω / R

    coriolis = BetaPlane(latitude = 0, β = β)

    # Vertical mixing: tiny background (preserve the thermocline) PLUS a
    # convective-adjustment safety valve.  A hydrostatic model cannot represent
    # convection, so when a column goes statically unstable (N² < 0) it must be
    # removed or the model blows up.  convective_κz fires ONLY where N² < 0, so
    # with a correct (restoring) surface heat flux it acts rarely instead of
    # homogenising the column.  convective_νz = 0 leaves the EUC momentum to the
    # biharmonic + Rayleigh drag.  Implicit, so it doesn't constrain Δt.
    vertical_mixing = ConvectiveAdjustmentVerticalDiffusivity(
        convective_κz = 1.0,
        convective_νz = 0.0,
        background_κz = 1e-5,
        background_νz = 1e-3)

    # Horizontal mixing: BIHARMONIC viscosity on MOMENTUM ONLY (κ defaults to 0).
    # Scale-selective — it damps the 2Δx grid-scale checkerboard noise (the
    # vertical stripes in w) by smoothing the velocity field that generates it,
    # while leaving the 1000s-of-km EUC almost untouched.
    #
    # NOTE: biharmonic is deliberately NOT applied to tracers.  Unlike Laplacian,
    # biharmonic diffusion is not sign-preserving and overshoots near sharp
    # gradients (the halocline); a κ here drove S negative → TEOS10 sqrt crash.
    # Tracers rely on WENO advection (overshoot-limited) + the small vertical κ.
    horizontal_mixing = HorizontalScalarBiharmonicDiffusivity(ν = 1e12)

    closure = (vertical_mixing, horizontal_mixing)


    model = HydrostaticFreeSurfaceModel(
        grid;
        free_surface        = ImplicitFreeSurface(),
        buoyancy            = buoyancy,
        tracers             = (:T, :S),
        coriolis            = coriolis,
        closure             = closure,
        momentum_advection  = WENO(order = 3),
        # High-order tracer advection: order-3 WENO is the most diffusive WENO and
        # was numerically mixing heat across the sharp thermocline (eroding the
        # eastern cold tongue's source water).  Order 7 has far less spurious
        # diapycnal mixing while staying non-oscillatory, so the thermocline
        # stays sharp.  (Needs grid halo ≥ 4 — set above.)
        tracer_advection    = WENO(order = 7),
        boundary_conditions = (u = u_bcs, v = v_bcs, T = T_bcs, S = S_bcs),
        forcing             = (u = u_forcing, v = v_forcing, T = T_forcing)
    )

    ################################################
    #### Initial conditions ####
    ################################################

    T_surface = 28.0
    T_mid     = 12.0
    T_deep    = 3.5

    z_thermocline    = -140.0
    δ_thermocline    = 55.0
    z_deep_transition = -900.0
    δ_deep_transition = 450.0

    Tᵢ(x, y, z) =
        T_deep +
        (T_mid - T_deep) *
        0.5 * (1 + tanh((z - z_deep_transition) / δ_deep_transition)) +
        (T_surface - T_mid) *
        0.5 * (1 + tanh((z - z_thermocline) / δ_thermocline))

    S_surface    = 34.6
    S_subsurface = 35.1
    S_deep       = 34.8

    z_halocline      = -120.0
    δ_halocline      = 70.0
    z_deep_freshening = -900.0
    δ_deep_freshening = 500.0

    Sᵢ(x, y, z) =
        S_deep +
        (S_subsurface - S_deep) *
        0.5 * (1 + tanh((z - z_deep_freshening) / δ_deep_freshening)) +
        (S_surface - S_subsurface) *
        0.5 * (1 + tanh((z - z_halocline) / δ_halocline))

    set!(model, u = 0, v = 0, w = 0, T = Tᵢ, S = Sᵢ)

    return model, grid, τx_parameters
end

################################################
#### Run model ####
################################################

function run_model(output_dir;
                   filename                 = "toy_model_test_equ_3d",
                   stop_time                = 5days,
                   plot_boundary_conditions = false,
                   plot_initial_conditions  = false)

    model, grid, τx_parameters = build_model()

    if plot_boundary_conditions
        plot_surface_wind_profile(grid, τx_parameters, output_dir)
    end

    if plot_initial_conditions
        plot_initial_conditions(model, grid, output_dir)
    end

    simulation = Simulation(model, Δt = 1, stop_time = stop_time)

    conjure_time_step_wizard!(simulation,
                              cfl    = 0.2,
                              max_Δt = 600)

    progress_message(sim) =
        @printf("Iteration: %04d, time: %s, Δt: %s, max(|w|) = %.1e m s⁻¹, wall time: %s\n",
                iteration(sim),
                prettytime(sim),
                prettytime(sim.Δt),
                maximum(abs, sim.model.velocities.w),
                prettytime(sim.run_wall_time))

    add_callback!(simulation, progress_message, IterationInterval(200))

    T    = model.tracers.T
    S    = model.tracers.S
    u, v, w = model.velocities
    dTdz = ∂z(T)

    save_interval = 48hours
    output_file   = joinpath(output_dir, filename * ".jld2")

    simulation.output_writers[:fields] =
        JLD2Writer(model,
                   (; T, S, u, v, w, dTdz),
                   filename           = output_file,
                   schedule           = TimeInterval(save_interval),
                   overwrite_existing = true)

    ################################################
    #### Momentum spin-down diagnostics ####
    ################################################
    # Volume-integrated zonal momentum split into the upper ocean (where the
    # EUC lives) and the deep ocean (below the Rayleigh-drag region).  U_deep
    # should relax toward ~0 and stay there if the deep ocean is spinning down;
    # if it keeps growing, increase λ_deep / raise z_damp_top.  U_deep_abs uses
    # |u| so opposite-sign deep flows can't cancel and hide residual momentum.
    z_diag = -300.0

    deep_mask_field  = Field{Center, Center, Center}(grid)
    upper_mask_field = Field{Center, Center, Center}(grid)
    set!(deep_mask_field,  (x, y, z) -> z <  z_diag ? 1.0 : 0.0)
    set!(upper_mask_field, (x, y, z) -> z >= z_diag ? 1.0 : 0.0)

    U_deep      = Field(Integral(u * deep_mask_field))        # ∫ u dV,   z < -300 m
    U_upper     = Field(Integral(u * upper_mask_field))       # ∫ u dV,   z ≥ -300 m
    U_deep_abs  = Field(Integral(abs(u) * deep_mask_field))   # ∫ |u| dV, z < -300 m
    KE_deep     = Field(Integral(u^2 * deep_mask_field))      # ∫ u² dV,  z < -300 m

    ts_interval = 24hours
    ts_file     = joinpath(output_dir, filename * "_momentum.jld2")

    simulation.output_writers[:momentum] =
        JLD2Writer(model,
                   (; U_deep, U_upper, U_deep_abs, KE_deep),
                   filename           = ts_file,
                   schedule           = TimeInterval(ts_interval),
                   overwrite_existing = true)

    run!(simulation)

    return model, grid, τx_parameters, output_file
end