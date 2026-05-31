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
    lat_edge = p.Ly / 2 / 111e3
    northern_peak = exp(-((lat - p.jet_lat) / p.jet_width)^2)
    southern_peak = exp(-((lat + p.jet_lat) / p.jet_width)^2)
    peak_shape = max(northern_peak, southern_peak)
    equator_shape = exp(-(lat / 8.0)^2)
    shape = peak_shape + 0.9 * equator_shape
    edge_taper = cos(π * lat / (2 * lat_edge))^2
    shape = shape * edge_taper
    shape = shape / 1.2
    return p.u₁₀_trades * shape
end

@inline function spinup_ramp(t, p)
    return min(1.0, t / p.spinup_time)
end

@inline function τx_profile(x, y, t, p)
    u10 = surface_u10_profile(y, p)
    τ = p.ρₐ / p.ρₒ * p.cᴰ * u10 * abs(u10)
    return τ * spinup_ramp(t, p)
end

################################################
#### North and south temperature fluxes ####
################################################

@inline function north_T_flux(x, z, p)
    depth = -z
    xnorm = x / p.Lx
    x_shape = 1.0 + 0.25 * sin(2π * xnorm) + 0.10 * sin(6π * xnorm + 0.7)
    surface_positive  =  1.8 * exp(-(depth / 20)^2)
    shallow_negative  = -1.3 * exp(-((depth - 80) / 55)^2)
    deeper_positive   =  0.4 * exp(-((depth - 100) / 150)^2)
    return p.surface_amp * x_shape * (surface_positive + shallow_negative + deeper_positive)
end

@inline function south_T_flux(x, z, p)
    depth = -z
    xnorm = x / p.Lx
    x_shape = 1.0 + 0.15 * sin(2π * xnorm + 0.5) + 0.05 * sin(6π * xnorm + 1.2)
    surface_negative        = -1.8 * exp(-(depth / 20)^2)
    shallow_positive        =  1.3 * exp(-((depth - 80) / 55)^2)
    shallow_negative        = -0.22 * exp(-((depth - 320) / 120)^2)
    shallow_positive_deeper =  0.28 * exp(-((depth - 430) / 70)^2)
    return p.surface_amp * x_shape * (surface_negative + shallow_positive + shallow_negative + shallow_positive_deeper)
end

@inline north_T_flux_bc(x, z, t, p) = north_T_flux(x, z, p)
@inline south_T_flux_bc(x, z, t, p) = south_T_flux(x, z, p)

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

@inline sponge_u(x, y, z, t, u, p) = -p.λ₀ * sponge_mask(x, y, p) * u
@inline sponge_v(x, y, z, t, v, p) = -p.λ₀ * sponge_mask(x, y, p) * v

################################################
#### Salinity boundary condition ####
################################################

@inline function Jˢ(x, y, t, S, p)
    evaporation_flux = -p.evaporation_rate * S
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
                           halo     = (3, 3, 3),
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

    Q   = 200
    cᴾ  = 3991
    Jᵀ  = Q / (ρₒ * cᴾ)

    north_T_flux_parameters = (;
        Lx,
        surface_amp = 5e-5
    )

    south_T_flux_parameters = (;
        Lx,
        surface_amp = 5e-5
    )

    T_bcs = FieldBoundaryConditions(
        top    = FluxBoundaryCondition(Jᵀ),
        bottom = FluxBoundaryCondition(0),
        north  = FluxBoundaryCondition(north_T_flux_bc, parameters = north_T_flux_parameters),
        south  = FluxBoundaryCondition(south_T_flux_bc, parameters = south_T_flux_parameters)
    )

    ################################################
    #### Salinity boundary conditions ####
    ################################################

    evaporation_rate  = 1e-8
    restore_rate      = 5e-7
    S_surface_target  = 34.6

    salinity_flux_parameters = (; evaporation_rate, restore_rate, S_surface_target)

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

    u_bcs = FieldBoundaryConditions(
        top    = FluxBoundaryCondition(τx_profile, parameters = τx_parameters),
        bottom = FluxBoundaryCondition(0)
    )

    v_bcs = FieldBoundaryConditions(
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

    sponge_params = (;
        sponge_width_x,
        sponge_width_y,
        Lx,
        Ly,
        λ₀,
        east_west_strength,
        north_south_strength
    )

    u_forcing = Forcing(sponge_u; field_dependencies = :u, parameters = sponge_params)
    v_forcing = Forcing(sponge_v; field_dependencies = :v, parameters = sponge_params)

    ################################################
    #### Model ####
    ################################################

    Ω = 7.2921e-5
    R = 6.371e6
    β = 2Ω / R

    coriolis = BetaPlane(latitude = 0, β = β)

    closure = ScalarDiffusivity(ν = 1e-4, κ = 1e-5)


    model = HydrostaticFreeSurfaceModel(
        grid;
        free_surface        = ImplicitFreeSurface(),
        buoyancy            = buoyancy,
        tracers             = (:T, :S),
        coriolis            = coriolis,
        closure             = closure,
        momentum_advection  = WENO(order = 3),
        tracer_advection    = WENO(order = 3),
        boundary_conditions = (u = u_bcs, v = v_bcs, T = T_bcs, S = S_bcs),
        forcing             = (u = u_forcing, v = v_forcing)
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

    return model, grid, τx_parameters, north_T_flux_parameters, south_T_flux_parameters
end

################################################
#### Run model ####
################################################

function run_model(output_dir;
                   filename                 = "toy_model_test_equ_3d",
                   stop_time                = 5days,
                   plot_boundary_conditions = false,
                   plot_initial_conditions  = false)

    model, grid, τx_parameters,
    north_T_flux_parameters,
    south_T_flux_parameters = build_model()

    if plot_boundary_conditions
        plot_surface_wind_profile(grid, τx_parameters, output_dir)
        plot_north_south_T_flux(north_T_flux, south_T_flux,
                                north_T_flux_parameters,
                                south_T_flux_parameters,
                                output_dir)
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

    run!(simulation)

    return model, grid, τx_parameters,
           north_T_flux_parameters,
           south_T_flux_parameters,
           output_file
end