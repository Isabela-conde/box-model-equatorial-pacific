using Oceananigans.Units: days
using Dates

include("toy_model.jl")
include("toy_plots.jl")

filename = "toy_model"

timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
output_dir = joinpath("runs", filename * "_" * timestamp)

mkpath(output_dir)
script_dir = @__DIR__

cp(joinpath(script_dir, "toy_model.jl"),
   joinpath(output_dir, "toy_model_cp.jl");
   force = true)

# cp(joinpath(script_dir, "toy_plots.jl"),
#    joinpath(output_dir, "toy_plots.jl");
#    force = true)


# cp(joinpath(script_dir, "run_toy_model.jl"),
#    joinpath(output_dir, "run_toy_model.jl");
#    force = true)


model,
grid,
τx_parameters,
north_T_flux_parameters,
south_T_flux_parameters,
filepath =
run_model(output_dir;
          filename = filename,
          stop_time = 1000days)




plot_surface_wind_profile(grid,
                         τx_parameters,
                         output_dir)

# plot_equator_full_depth(filepath,
#                        grid,
#                        output_dir)

# plot_equator_top_400m(filepath,
#                      grid,
#                      output_dir)

# output_dir = "/Users/isabelaconde/Documents/GitHub/toy_model_surface_wind_flux_increase/runs/toy_model_2026-05-15_115819/"
# filepath = joinpath(output_dir, filename * ".jld2")

# u = FieldTimeSeries(filepath, "u")
# grid = u.grid


plot_ymean_full_depth(filepath,
                     grid,
                     output_dir)

plot_ymean_top_400m(filepath,
                   grid,
                   output_dir)

plot_surface_video(filepath,
                  grid,
                  output_dir)

plot_equatorial_band_cross_section(filepath,
                                  grid,
                                  output_dir)


# Ly = grid.Ly

# τx_parameters = (;
#     u₁₀_trades = -10.0,
#     jet_width = 2.5,
#     jet_lat = 6.0,
#     Ly = grid.Ly,
#     ρₐ = 1.225,
#     ρₒ = 1025.0,
#     cᴰ = 1.3e-3,
# )


# plot_surface_wind_boundary(grid,
#                            τx_parameters,
#                            output_dir)

println("Finished run.")
println("Outputs saved in: ", output_dir)

