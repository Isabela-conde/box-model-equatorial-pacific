using Oceananigans
using Oceananigans.Units: hours
using CairoMakie
using Printf
using Statistics


Nx = 140
Ny = 40
Nz = 60

Lx = 1.5e7
Ly = 2.22e6
Lz = 2000

refinement = 3.5
stretching = 18

h(k) = (k - 1) / Nz

ζ₀(k) = 1 + (h(k) - 1) / refinement

Σ(k) = (1 - exp(-stretching * h(k))) /
       (1 - exp(-stretching))

z_interfaces(k) = Lz * (ζ₀(k) * Σ(k) - 1)


# Generating function
z_interfaces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

grid = RectilinearGrid(CPU(),
                       size = (Nx, Nx, Nz),
                       x = (0, Lx),
                       y = (0, Ly),
                       z = z_interfaces)




# plotting 

Δz = Array(zspacings(grid, Center()))[1, 1, :]
zc = Array(znodes(grid, Center()))

fig = Figure(size = (1200, 800))
ax = Axis(fig[1, 1], xlabel = "Vertical spacing (m)", ylabel = "z (m)")

lines!(ax, Δz, zc)
scatter!(ax, Δz, zc)

fig


####################################### temp advec for north and south fluxes

Lx = 1.5e7
Lz = 2000
surface_amp = 1e-8

@inline function north_T_flux(x, z, p)
    depth = -z
    xnorm = x / p.Lx

    x_shape =
        1.0 +
        0.25 * sin(2π * xnorm) +
        0.10 * sin(6π * xnorm + 0.7)

    surface_positive =
        1.8 * exp(-(depth / 20)^2)

    shallow_negative =
        -1.3 * exp(-((depth - 80) / 55)^2)

    deeper_positive =
        0.4 * exp(-((depth - 100) / 150)^2)

    return p.surface_amp * x_shape *
           (surface_positive + shallow_negative + deeper_positive)
end

@inline function south_T_flux(x, z, p)
    depth = -z
    xnorm = x / p.Lx

    x_shape =
        1.0 +
        0.15 * sin(2π * xnorm + 0.5) +
        0.05 * sin(6π * xnorm + 1.2)

    surface_negative =
        -1.8 * exp(-(depth / 20)^2)

    # right_half_mask =
    #     exp(-((xnorm - 0.72) / 0.22)^2)

    # middle_positive =
    #     1.8 *
    #     right_half_mask *
    #     exp(-((depth - 100) / 85)^2)

    shallow_positive =
        1.3 * exp(-((depth - 80) / 55)^2)

    shallow_negative =
        -0.22 * exp(-((depth - 320) / 120)^2)

    shallow_positive_deeper =
        0.28 * exp(-((depth - 430) / 70)^2)

    vertical_structure =
        surface_negative +
        # middle_positive +
        shallow_positive +
        shallow_negative +
        shallow_positive_deeper

    return p.surface_amp * x_shape * vertical_structure
end

params = (; Lx = Lx, surface_amp = surface_amp)

xs = range(0, Lx, length = 300)
zs = range(-Lz, 0, length = 200)

north_field = [north_T_flux(x, z, params) for x in xs, z in zs]
south_field = [south_T_flux(x, z, params) for x in xs, z in zs]

clim = maximum(abs, vcat(vec(north_field), vec(south_field)))

fig = Figure(size = (1100, 450))

ax1 = Axis(fig[1, 1],
    xlabel = "x distance (km)",
    ylabel = "Depth (m)",
    title = "North T flux")

hm1 = heatmap!(ax1, xs ./ 1e3, -zs, north_field,
    colormap = Reverse(:RdBu),
    colorrange = (-clim, clim))

ax1.yreversed = true

ax2 = Axis(fig[1, 2],
    xlabel = "x distance (km)",
    ylabel = "Depth (m)",
    title = "South T flux")

hm2 = heatmap!(ax2, xs ./ 1e3, -zs, south_field,
    colormap = Reverse(:RdBu),
    colorrange = (-clim, clim))

ax2.yreversed = true

Colorbar(fig[1, 3], hm2, label = "T flux")

fig

# save("north_south_T_flux.png", fig)





######################################################


Ly = 2.22e6

u₁₀_equator = -2
u₁₀_trades  = -10
jet_lat     = 5.0
jet_width   = 2.5

cᴰ = 2e-3
ρₐ = 1.2
ρₒ = 1026

τx_parameters = (;
    u₁₀_equator,
    u₁₀_trades,
    jet_lat,
    jet_width,
    cᴰ,
    ρₐ,
    ρₒ,
    Ly
)

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

@inline function τx_profile(x, y, t, p)
    u10 = surface_u10_profile(y, p)
    return p.ρₐ / p.ρₒ * p.cᴰ * u10 * abs(u10)
end

ys = range(-Ly/2, Ly/2, length = 500)

lats = ys ./ 111e3

u10 = [
    surface_u10_profile(y, τx_parameters)
    for y in ys
]

fig = Figure(size = (700, 400))

ax = Axis(
    fig[1,1],
    xlabel = "Latitude",
    ylabel = "u₁₀ (m s⁻¹)",
    title = "Idealised zonal wind profile"
)

lines!(ax, lats, u10, linewidth = 4)

hlines!(ax, [0], linestyle = :dash)

fig

# save("surface_wind_profile.png", fig)



###################### salinity 


Ly = 2.22e6
flux_amp = 5e-7

salinity_flux_parameters = (; flux_amp, Ly)

@inline function Jˢ(x, y, t, S, p)

    y_norm = y / (p.Ly / 2)

    freshwater_pattern =
        -cos(π * y_norm)

    return freshwater_pattern * p.flux_amp * S
end

ys = range(-Ly/2, Ly/2, length = 500)
lats = ys ./ 111e3

S0 = 35.0

salt_flux = [
    Jˢ(0, y, 0, S0, salinity_flux_parameters)
    for y in ys
]

fig = Figure(size = (700, 400))

ax = Axis(
    fig[1, 1],
    xlabel = "Latitude",
    ylabel = "Surface salinity flux",
    title = "Surface salinity boundary condition"
)

lines!(ax, lats, salt_flux, linewidth = 4)
hlines!(ax, [0], linestyle = :dash)

fig

# save("surface_salinity_flux.png", fig)



