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


# relative_wind = true  → current-feedback (relative-wind) stress
# relative_wind = false → resting-ocean stress
# (set to false for now to isolate whether the relative wind is killing the SEC)
model,
grid,
τx_parameters,
filepath =
run_model(output_dir;
          filename = filename,
          stop_time = 400days,   # MUST run this long: the EUC/SEC goes through a
                                 # spin-up OVERSHOOT (~day 100-200 the EUC shoals to
                                 # ~-25m & the surface nearly flips eastward) and only
                                 # recovers to a strong westward SEC by ~day 250-300.
                                 # checkpoint_2 (the good run) did exactly this — short
                                 # 150-day runs stop mid-overshoot and look broken.
          relative_wind = false,
          poleward_sink = false)





          plot_surface_wind_profile(grid,
                         τx_parameters,
                         output_dir)

momentum_file = joinpath(output_dir, filename * "_momentum.jld2")
plot_momentum_diagnostics(momentum_file, output_dir)

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

################################################
#### Last-100-day analysis plots ####
################################################

# ȳ-mean cross sections (T, b, ∂b/∂z, u) with thermocline overlay:
plot_ymean_last100(filepath, output_dir; window_days = 10, depth_max = 700.0)

# Equatorial-band (|y| ≤ 250 km) cross sections:
plot_ymean_last100(filepath, output_dir; window_days = 10, depth_max = 700.0, y_band = 250e3)

# Meridional sections (latitude on the horizontal axis) at two classic EUC
# longitudes — 140°W and 110°W (idealised basin assumed to span 120°E→80°W):
plot_meridional_last100(filepath, output_dir; window_days = 10, depth_max = 700.0,
                        xloc = pacific_lon_to_x(-140.0), label = "≈140°W", tag = "140W")
plot_meridional_last100(filepath, output_dir; window_days = 10, depth_max = 700.0,
                        xloc = pacific_lon_to_x(-110.0), label = "≈110°W", tag = "110W")

# Surface fields (last-100-day mean):
plot_surface_last100(filepath, output_dir; window_days = 100)

# Surface heat budget — uses the function's default flux params, which track the
# current model (λ_T=7e-6, evaporation, 70 W/m² solar).  If you change the model's
# surface-flux params, update plot_surface_heat_budget's defaults (or pass them here).
plot_surface_heat_budget(filepath, output_dir)



println("Finished run.")
println("Outputs saved in: ", output_dir)

