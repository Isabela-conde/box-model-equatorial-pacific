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

    # Smooth off-equatorial trade-wind bump (no equatorial cusp).  The factor
    # (lat/jet_lat)² ·exp(1-(lat/jet_lat)²) is an even, C∞ function of latitude
    # that is exactly 0 at the equator (→ u₁₀_equator), exactly 1 at the trade
    # latitudes ±jet_lat (→ u₁₀_trades).
    r           = lat / p.jet_lat
    trade_shape = r^2 * exp(1 - r^2)

    # Edge taper: subtract a smooth, poleward-increasing term so the easterlies
    # WEAKEN toward the N/S walls but stay NON-ZERO there (unlike the checkpoint_2
    # cos²-to-zero taper).  (lat/lat_edge)⁴ ≈0 through the interior and = 1 at the
    # walls, where u₁₀_wall_drop is removed.  C∞ everywhere — no kink.
    lat_edge   = p.Ly / 2 / 111e3
    edge_taper = (lat / lat_edge)^4

    return p.u₁₀_equator + (p.u₁₀_trades - p.u₁₀_equator) * trade_shape -
           p.u₁₀_wall_drop * edge_taper
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

# NOTE: the N/S walls are INSULATING (zero T flux), consistent with the salinity walls
# and the bottom.  The checkpoint_2 wall heat-flux dipole was carried along until
# 2026-10-21 but contributed +0.0003 PW — nothing measurable — and is gone; the poleward
# heat sink (§ ns_band / sink_T) provides the actual net heat export.

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

# Surface heat flux = SST restoring ONLY (reverted to checkpoint_2's pure-restoring
# form).  Oceananigans convention: +ve flux = heat OUT, so this is a negative
# feedback toward the zonally-structured target T_target_profile(x) (warm pool W →
# cold tongue E).  The latent-cooling and penetrating-solar terms added after
# checkpoint_2 have been removed.
@inline function surface_T_flux_bc(x, y, t, T, p)
    return p.λ_T * (T - T_target_profile(x, p))
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

# DISCRETE-FORM versions — the ones the model uses (same physics as above).  Like the
# sink (see sink_T_discrete), these index model_fields by NAME so they cannot receive
# the wrong field through Oceananigans 0.96.35's continuous `field_dependencies`
# forcing bug.  u lives on x-faces and v on y-faces, so the node functions use the
# matching locations.
@inline function sponge_u_discrete(i, j, k, grid, clock, model_fields, p)
    x = Oceananigans.Grids.xnode(i, j, k, grid, Face(), Center(), Center())
    y = Oceananigans.Grids.ynode(i, j, k, grid, Face(), Center(), Center())
    z = Oceananigans.Grids.znode(i, j, k, grid, Face(), Center(), Center())
    return -(p.λ₀ * sponge_mask(x, y, p) + p.λ_deep * deep_mask(z, p)) * @inbounds(model_fields.u[i, j, k])
end
@inline function sponge_v_discrete(i, j, k, grid, clock, model_fields, p)
    x = Oceananigans.Grids.xnode(i, j, k, grid, Center(), Face(), Center())
    y = Oceananigans.Grids.ynode(i, j, k, grid, Center(), Face(), Center())
    z = Oceananigans.Grids.znode(i, j, k, grid, Center(), Face(), Center())
    return -(p.λ₀ * sponge_mask(x, y, p) + p.λ_deep * deep_mask(z, p)) * @inbounds(model_fields.v[i, j, k])
end

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

# Poleward heat sink (the only T body forcing): relax T toward a cool extratropical
# reference in the N/S bands.  Solar penetration was removed with the heat-flux
# revert.  Attached only when poleward_sink=true (currently off).
# Depth mask for the sink: ≈1 above ~150 m, 0.5 at z_sink_mid, ≈0 below ~450 m.
# Without it the sink chilled the 300–900 m band water (−6.6 PW of the −8.5 PW total
# in run 2026-10-20_150728), which is the EUC source water → 14 °C cold tongue and a
# basin cooling at −2.7 PW.  Confining the sink to the upper ocean leaves the
# poleward export of the cold-tongue heat uptake, which is its job.
@inline sink_depth_mask(z, p) = 0.5 * (1 + tanh((z - p.z_sink_mid) / p.d_sink))

@inline sink_T(x, y, z, t, T, p) =
    -p.λ_T_sink * ns_band_mask(y, p) * sink_depth_mask(z, p) * (T - T_extratropical_ref(z, p))

# DISCRETE-FORM version of the same sink — this is the one the model uses.
# ⚠ Oceananigans 0.96.35 (Gadi) hands a continuous `Forcing(sink_T; field_dependencies = :T)`
# the WRONG field (≈0, apparently η) instead of T, so the sink evaluated
# −λ·(0 − T_ref) = +λ·T_ref: a constant ~+15 PW HEATER (600-day run 2026-10-20_114944:
# ocean warmed 250 W/m², 36 °C wall strips, thermocline sinking).  The discrete form
# indexes model_fields.T by NAME, so it is immune to that bug; same physics as sink_T.
# (0.107 on the laptop evaluates the continuous form correctly; both forms agree there.)
@inline function sink_T_discrete(i, j, k, grid, clock, model_fields, p)
    # ynode/znode (singular) are not exported at top level → qualify (works in 0.96 and 0.107)
    y = Oceananigans.Grids.ynode(i, j, k, grid, Center(), Center(), Center())
    z = Oceananigans.Grids.znode(i, j, k, grid, Center(), Center(), Center())
    return -p.λ_T_sink * ns_band_mask(y, p) * sink_depth_mask(z, p) *
           (@inbounds(model_fields.T[i, j, k]) - T_extratropical_ref(z, p))
end

################################################
#### Salinity boundary condition ####
################################################

# Same convention as the T flux: a POSITIVE top flux REMOVES salt.  So the restoring
# must be +restore_rate·(S − S_target): S above target → positive flux → salt out →
# negative feedback.  ⚠ Until 2026-10-21 this had a minus sign (anti-restoring, the
# same error the T flux once had): basin-mean SSS drifted 34.61 → 35.26 over the
# 5-year spin-up without ever approaching 34.6.  With the sign fixed, evaporation and
# restoring balance at SSS ≈ 35.3, so the fix stops the drift without a shock.
@inline function Jˢ(x, y, t, S, p)
    evaporation_flux = -p.evaporation_rate * S   # uniform evaporation (checkpoint_2); -ve → S increases
    restoring_flux   = +p.restore_rate * (S - p.S_surface_target)
    return evaporation_flux + restoring_flux
end

################################################
#### Build model ####
################################################

function build_model(; relative_wind = true, poleward_sink = true,
                       κ₀_Ri = 1e-2,   # Ri-scheme max tracer diffusivity (branch experiments)
                       ν₀_Ri = 1e-2,   # Ri-scheme max shear viscosity    (branch experiments)
                       cᴰ    = 2e-3)   # wind drag coefficient; 1.2e-3 is the standard value (branch experiments)

    Nx = 600         # Δx ≈ 25 km (was 140 → 107 km)
    Ny = 150         # ±15° at Δy ≈ 22 km (was 60 → 55 km)
    Nz = 60

    Lx = 1.5e7
    Ly = 3.33e6      # ±15° (was 2.22e6 = ±10°); 500 km sponges/sink bands now sit at ~10.5–15°
    Lz = 2000

    refinement = 3.5
    stretching = 18

    h(k)  = (k - 1) / Nz
    ζ₀(k) = 1 + (h(k) - 1) / refinement
    Σ(k)  = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))
    z_interfaces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

    grid = RectilinearGrid(GPU(),
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
    # λ_T = 2e-5 ≈ 80 W m⁻² K⁻¹ → restoring timescale over a 50 m layer ≈ 29 days.
    # Restored to checkpoint_2's strong value: weak λ_T let the cold tongue crash
    # (→ eastern convection) and the warm pool drift (→ flat thermocline).  Target
    # SST is zonally structured: warm pool (T_warm) west, cold tongue (T_cold) east.
    # 2026-10-20: 2e-5 → 4e-5 (≈160 W m⁻² K⁻¹).  At 2e-5 the cold tongue sat at 16 °C
    # vs the 24 °C target (eastern-wall upwelling out-pumped the restoring; 800-day run
    # 2026-10-20_182149).  Doubling pins SST harder without touching the dynamics; the
    # price is a damped SST response in the warming experiments.
    λ_T         = 4e-5
    T_warm      = 29.5     # western warm-pool target
    T_cold      = 24.0     # eastern cold-tongue target
    x_front     = 0.6      # front location as fraction of basin width (W→E)
    front_width = 0.18     # front sharpness (fraction of Lx)

    surface_T_flux_parameters = (; λ_T, T_warm, T_cold, x_front, front_width, Lx)

    surface_T_bc = FluxBoundaryCondition(surface_T_flux_bc,
                                         field_dependencies = :T,
                                         parameters         = surface_T_flux_parameters)

    # N/S walls insulating (default no-flux).  The checkpoint_2 wall heat-flux dipole
    # was removed 2026-10-21: its net was +0.0003 PW and its per-cell tendency ~4000×
    # smaller than the SST restoring, i.e. it did nothing measurable.
    T_bcs = FieldBoundaryConditions(
        top    = surface_T_bc,
        bottom = FluxBoundaryCondition(0)
    )

    ################################################
    #### Salinity boundary conditions ####
    ################################################

    evaporation_rate = 1e-8     # uniform evaporation (checkpoint_2)
    restore_rate     = 5e-7
    S_surface_target = 34.6

    salinity_flux_parameters = (; evaporation_rate, restore_rate, S_surface_target)

    evaporation_bc = FluxBoundaryCondition(Jˢ,
                                           field_dependencies = :S,
                                           parameters         = salinity_flux_parameters)

    S_bcs = FieldBoundaryConditions(top = evaporation_bc)

    ################################################
    #### Wind boundary condition ####
    ################################################

    # Smooth-bump + wall-taper wind (non-zero at the walls).  Equator dip
    # u₁₀_equator, off-equatorial trade peak u₁₀_trades at ±jet_lat, easterlies
    # tapering to below the equatorial value (but NOT zero) at the N/S walls.
    u₁₀_equator   = 7.0    # easterly wind speed at the equator (dip)
    u₁₀_trades    = 8.0    # easterly maximum at the trade latitudes ±jet_lat
    jet_lat       = 10   # trade peak latitude (±10°, retuned for the ±15° domain)
    u₁₀_wall_drop = 4.0    # easterly reduction at the N/S walls (→ wall ≈ 3.6 m/s at ±15°)
    # cᴰ is a build_model keyword (default 2e-3).  With u₁₀ = 7 m/s that gives
    # τ_eq = 0.118 N m⁻², ~2× the observed annual mean at 140W; cᴰ = 1.2e-3 → 0.071.
    ρₐ          = 1.2
    spinup_time = 10days

    τx_parameters = (;
        u₁₀_equator,
        u₁₀_trades,
        jet_lat,
        u₁₀_wall_drop,
        cᴰ,
        ρₐ,
        ρₒ,
        Ly,
        spinup_time
    )

    # Wind stress.  Toggle with the `relative_wind` keyword to build_model:
    #   relative_wind = true  → current-feedback stress (τ depends on surface u,v;
    #                           τy opposes the ocean's own v — eddy damping)
    #   relative_wind = false → resting-ocean stress (τx_profile only; v has no
    #                           top stress).  At rest the two are identical.
    if relative_wind
        u_bcs = FieldBoundaryConditions(
            top    = FluxBoundaryCondition(τx_relative, field_dependencies = (:u, :v),
                                           parameters = τx_parameters),
            bottom = FluxBoundaryCondition(0))
        v_bcs = FieldBoundaryConditions(
            top    = FluxBoundaryCondition(τy_relative, field_dependencies = (:u, :v),
                                           parameters = τx_parameters),
            bottom = FluxBoundaryCondition(0))
    else
        u_bcs = FieldBoundaryConditions(
            top    = FluxBoundaryCondition(τx_profile, parameters = τx_parameters),
            bottom = FluxBoundaryCondition(0))
        v_bcs = FieldBoundaryConditions(
            bottom = FluxBoundaryCondition(0))
    end

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

    # discrete_form: see sponge_u_discrete — avoids the 0.96.35 field_dependencies bug.
    u_forcing = Forcing(sponge_u_discrete; discrete_form = true, parameters = sponge_params)
    v_forcing = Forcing(sponge_v_discrete; discrete_form = true, parameters = sponge_params)

    # Poleward heat sink: relax T toward a cool extratropical profile in the N/S
    # bands so the net tropical heat uptake is exported and the basin equilibrates.
    λ_T_sink   = 1 / 150days
    T_surf_ref = 22.0      # cool extratropical surface (°C)
    T_deep_ref = 3.5       # deep reference (matches initial deep T → no deep nudging)
    z_ref      = -150.0    # reference thermocline depth (m)
    d_ref      = 120.0     # reference thermocline thickness (m)
    z_sink_mid = -300.0    # sink depth mask = 0.5 here (see sink_depth_mask)
    d_sink     =   75.0    # taper width (m): mask ≈0.98 at 150 m, ≈0.02 at 450 m

    sink_params = (;
        sponge_width_y,
        Ly,
        λ_T_sink,
        T_surf_ref,
        T_deep_ref,
        z_ref,
        d_ref,
        z_sink_mid,
        d_sink
    )

    # Toggle the poleward heat sink with the `poleward_sink` keyword to build_model.
    # discrete_form: see sink_T_discrete — the continuous `field_dependencies = :T`
    # form gets the wrong field in Oceananigans 0.96.35 and turns the sink into a heater.
    T_forcing = Forcing(sink_T_discrete; discrete_form = true, parameters = sink_params)

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

    # Richardson-number-dependent VERTICAL mixing (Oceananigans built-in), ADDED on
    # top of the constant background above.  It raises κ and ν where the gradient
    # Richardson number Ri = N²/(u_z²+v_z²) is low (strong shear / weak stratification
    # — the surface layer and EUC flanks) and tapers them to ~0 in the sharply
    # stratified thermocline, so mixing is concentrated where shear instability
    # actually mixes.  Form: κ = κ₀·τ(Ri), ν = ν₀·τ(Ri) with the tanh taper
    # τ = ½(1 - tanh((Ri-Ri₀)/Riᵟ)).  κᶜᵃ = 0 here because the convective valve is
    # already provided by the ConvectiveAdjustmentVerticalDiffusivity above (avoid
    # double-counting N²<0).  Vertically implicit → does NOT constrain Δt.
    # ⚠ Oceananigans flags this closure as experimental / parameters uncalibrated;
    # ν₀, κ₀, Ri₀, Riᵟ are starting values to TUNE against the EUC/thermocline.
    ri_mixing = RiBasedVerticalDiffusivity(; warning = false,
                                             ν₀  = ν₀_Ri,  # max shear viscosity (low Ri); build_model keyword
                                             # max shear diffusivity (low Ri).  Now a build_model
                                             # keyword: with κ₀ = ν₀ (Pr = 1) the scheme mixed the
                                             # restoring heat 75 m down at the equator → N² ~3e-5,
                                             # Ri 0.1–0.25 → ν ~6e-3 → EUC momentum reached the
                                             # surface, no SEC (spin-up 2026-10-21_005204).  PP81
                                             # has κ ≈ ν/(1+5Ri) ≈ ν/3 here → test κ₀ = 3e-3.
                                             κ₀  = κ₀_Ri,
                                             κᶜᵃ = 0.0,    # convection handled by CAVD above
                                             Ri₀ = 0.1,    # taper centre (Richardson number)
                                             Riᵟ = 0.4,    # taper width
                                             # Smooth Ri horizontally BEFORE computing κ,ν: breaks the
                                             # noisy-shear → noisy-Ri → noisy-κ feedback that made the
                                             # 2Δx equatorial lines at coarse Δx.  Not exported → qualify.
                                             horizontal_Ri_filter = Oceananigans.TurbulenceClosures.FivePointHorizontalFilter())

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

    closure = (vertical_mixing, ri_mixing, horizontal_mixing)
    # To revert to the pre-Ri (checkpoint_2) closure, drop ri_mixing:
    # closure = (vertical_mixing, horizontal_mixing)


    model = HydrostaticFreeSurfaceModel(;
        grid,   # keyword form: the only one Oceananigans 0.96 accepts (newer versions also allow positional)
        free_surface        = ImplicitFreeSurface(),
        buoyancy            = buoyancy,
        tracers             = (:T, :S),
        coriolis            = coriolis,
        closure             = closure,
        momentum_advection  = WENO(order = 3),
        # Tracer advection: WENO order 5 (lowered from 7).  Order 7 is the least
        # dissipative WENO and left the equatorial cold-tongue front (advected by the
        # SEC at dx=107 km) ringing with undamped 2Δx noise (surface T zig-zagging
        # ±1-2.5 °C, corrupting the shallow eastern thermocline).  Order 5 adds more
        # built-in dissipation to damp that advective dispersion, while staying far
        # less diffusive than order 3 (which originally eroded the thermocline).
        # Watch: (a) does the eastern zig-zag collapse, (b) does the deep western
        # thermocline stay sharp.  (Needs grid halo ≥ 3; halo is 4.)
        tracer_advection    = WENO(order = 5),
        boundary_conditions = (u = u_bcs, v = v_bcs, T = T_bcs, S = S_bcs),
        forcing             = poleward_sink ?
                              (u = u_forcing, v = v_forcing, T = T_forcing) :
                              (u = u_forcing, v = v_forcing)
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
                   relative_wind            = true,
                   poleward_sink            = true,
                   plot_boundary_conditions = false,
                   plot_initial_conditions  = false,
                   save_interval            = 30days,   # 3-D time-mean output cadence
                   pickup                   = false,    # false = from rest; or path to a checkpoint file
                   κ₀_Ri                    = 1e-2,     # Ri-scheme tracer diffusivity (see build_model)
                   ν₀_Ri                    = 1e-2,     # Ri-scheme shear viscosity
                   cᴰ                       = 2e-3)     # wind drag coefficient

    model, grid, τx_parameters = build_model(; relative_wind = relative_wind,
                                               poleward_sink = poleward_sink,
                                               κ₀_Ri = κ₀_Ri, ν₀_Ri = ν₀_Ri, cᴰ = cᴰ)

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

    # ---- 3-D fields: TIME-AVERAGED and sparse (the big file) ----------------------
    # Instantaneous 3-D snapshots every 48 h were ~259 MB each at 600×150×60 (6 fields)
    # → ~78 GB for a 600-day run.  Instead write 10-day MEANS every 10 days: 60 snapshots
    # for 600 d (~15 GB).  Means are what the equilibrium / thermocline / EUC diagnostics
    # want anyway.  Same fields and names, so every existing plot still works unchanged.
    # (The eastern 2Δx noise is a STANDING pattern, so it survives averaging — the
    # line-check still works on this output.)  save_interval is a run_model keyword:
    # 10days for ≤800-day runs, 30days for multi-year spin-ups (~10 GB per 5 years).
    output_file   = joinpath(output_dir, filename * ".jld2")

    simulation.output_writers[:fields] =
        JLD2Writer(model,
                   (; T, S, u, v, w, dTdz),
                   filename           = output_file,
                   schedule           = AveragedTimeInterval(save_interval; window = save_interval),
                   overwrite_existing = true)

    # ---- 2-D slices: small, INSTANTANEOUS, frequent (every 2 days) ----------------
    # These carry the high-time-resolution diagnostics so the 3-D file doesn't have to.
    # All four together are ~1.5 GB per 600 d.
    slice_interval = 2days
    Nx, Ny, Nz     = size(grid)

    # surface layer (top cell) — feeds the surface video / surface maps.
    # dTdz included so plot_surface_video's ∂T/∂z panel works off this file.
    simulation.output_writers[:surface] =
        JLD2Writer(model,
                   (; T, S, u, v, w, dTdz),
                   filename           = joinpath(output_dir, filename * "_surface.jld2"),
                   indices            = (:, :, Nz),
                   schedule           = TimeInterval(slice_interval),
                   overwrite_existing = true)

    # meridional (y–z) sections at ≈140°W and ≈110°W.  The idealised basin spans
    # 120°E→80°W, so 140°W ≈ 0.625·Lx and 110°W ≈ 0.81·Lx (same mapping as
    # pacific_lon_to_x in toy_plots.jl).
    i140 = round(Int, 0.625 * Nx)
    i110 = round(Int, 0.81  * Nx)

    simulation.output_writers[:sec140W] =
        JLD2Writer(model,
                   (; T, S, u, v, w),
                   filename           = joinpath(output_dir, filename * "_140W.jld2"),
                   indices            = (i140, :, :),
                   schedule           = TimeInterval(slice_interval),
                   overwrite_existing = true)

    simulation.output_writers[:sec110W] =
        JLD2Writer(model,
                   (; T, S, u, v, w),
                   filename           = joinpath(output_dir, filename * "_110W.jld2"),
                   indices            = (i110, :, :),
                   schedule           = TimeInterval(slice_interval),
                   overwrite_existing = true)

    # zonal (x–z) y-MEAN section — thermocline tilt and the EUC along the basin.
    # (Equatorial-BAND means are still computed by the plots from the 3-D file.)
    T_ymean = Field(Average(T; dims = 2))
    u_ymean = Field(Average(u; dims = 2))

    simulation.output_writers[:ymean] =
        JLD2Writer(model,
                   (; T = T_ymean, u = u_ymean),
                   filename           = joinpath(output_dir, filename * "_ymean.jld2"),
                   schedule           = TimeInterval(slice_interval),
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

    # ---- Checkpointer: full model state every 100 days, latest kept ----------------
    # Enables spin-up → branch experiments: a later run passes this file as `pickup`
    # and continues from this state.  Saves u, v, T, S, η, the AB2 tendencies and the
    # clock (~0.5 GB per file; cleanup = true deletes the previous one).
    simulation.output_writers[:checkpointer] =
        Checkpointer(model;
                     schedule = TimeInterval(100days),
                     dir      = output_dir,
                     prefix   = "checkpoint",
                     cleanup  = true)

    # pickup = false → start from rest (initial conditions above).
    # pickup = "<run>/checkpoint_iteration<N>.jld2" → restore that state AND clock, so
    # stop_time is then ABSOLUTE (spin-up end + branch length).  Forcing/BC parameters
    # may differ from the checkpointed run (that is the perturbation); grid and closures
    # must not.  The file path must be explicit: pickup = true would look in THIS run's
    # (empty) output_dir.
    run!(simulation; pickup)

    return model, grid, τx_parameters, output_file
end