using Oceananigans.Units: days
using Dates
using CUDA   # GPU run (per the Gadi guide); needs Pkg.add("CUDA") in this project env

include("toy_model.jl")
include("toy_plots.jl")

filename = "toy_model"

# ---------------------------------------------------------------------------
# WHICH RUN?  Chosen at submit time so several branches can sit in the queue
# together without editing this file between submissions:
#     qsub -v BRANCH=control   run_toy_model.sh
#     qsub -v BRANCH=nu5e-3    run_toy_model.sh
#     qsub run_toy_model.sh                      (no BRANCH → "spinup")
#   spinup     : 5-year run from rest (1825 d, 30-day means, checkpoint every 100 d)
#   control    : 1-year branch from the spin-up checkpoint, nothing changed
#   kappa3e-3  : branch, Ri-scheme tracer diffusivity κ₀ 1e-2 → 3e-3  (done: no SEC change)
#   nu5e-3     : branch, Ri-scheme shear viscosity    ν₀ 1e-2 → 5e-3  (SEC test 2)
#   cd1.2e-3   : branch, wind drag cᴰ 2e-3 → 1.2e-3 (τ_eq 0.118 → 0.071 N/m²) (SEC test 3)
# Branch stop_time is ABSOLUTE = checkpoint time + branch length.  The checkpointer
# fires every 100 days, so the 1825-day spin-up's last checkpoint is at day 1800.
# ---------------------------------------------------------------------------
branch      = get(ENV, "BRANCH", isempty(ARGS) ? "spinup" : ARGS[1])
spinup_ckpt = "runs/toy_model_2026-10-21_005204/checkpoint_iteration259762.jld2"   # model day 1800

# defaults = the spin-up configuration; each branch overrides what it tests
κ₀_Ri = 1e-2;  ν₀_Ri = 1e-2;  cᴰ = 2e-3
if branch == "spinup"
    pickup = false;        stop_time = 1825days;             save_interval = 30days
elseif branch == "control"
    pickup = spinup_ckpt;  stop_time = 1800days + 365days;   save_interval = 10days
elseif branch == "kappa3e-3"
    pickup = spinup_ckpt;  stop_time = 1800days + 365days;   save_interval = 10days;  κ₀_Ri = 3e-3
elseif branch == "nu5e-3"
    pickup = spinup_ckpt;  stop_time = 1800days + 365days;   save_interval = 10days;  ν₀_Ri = 5e-3
elseif branch == "cd1.2e-3"
    pickup = spinup_ckpt;  stop_time = 1800days + 365days;   save_interval = 10days;  cᴰ = 1.2e-3
else
    error("unknown BRANCH = \"$branch\" (use spinup, control, kappa3e-3, nu5e-3 or cd1.2e-3)")
end
println("BRANCH = $branch | pickup = $pickup | stop_time = $(stop_time/days) days | κ₀_Ri = $κ₀_Ri | ν₀_Ri = $ν₀_Ri | cᴰ = $cᴰ")

timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
output_dir = joinpath("runs", filename * "_" * branch * "_" * timestamp)

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
# pickup / stop_time / save_interval / κ₀_Ri come from the BRANCH block above.
model,
grid,
τx_parameters,
filepath =
run_model(output_dir;
          filename = filename,
          stop_time = stop_time,
          save_interval = save_interval,
          pickup = pickup,
          κ₀_Ri = κ₀_Ri,
          ν₀_Ri = ν₀_Ri,
          cᴰ = cᴰ,
          relative_wind = false,
          poleward_sink = true)


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

# Surface video now reads the 2-day INSTANTANEOUS surface-slice file written by
# run_model (the 3-D file `filepath` is 10-day means, so it has far fewer frames).
surface_file = joinpath(output_dir, filename * "_surface.jld2")
plot_surface_video(surface_file,
                  grid,
                  output_dir)

plot_equatorial_band_cross_section(filepath,
                                  grid,
                                  output_dir)

################################################
#### Last-100-day analysis plots ####
################################################

# ȳ-mean cross sections (T, b, ∂b/∂z, u) with thermocline overlay:
plot_ymean_last100(filepath, output_dir; window_days = 100, depth_max = 400.0)

# Equatorial-band (|y| ≤ 250 km) cross sections:
plot_ymean_last100(filepath, output_dir; window_days = 100, depth_max = 400.0, y_band = 250e3)

# Meridional sections (latitude on the horizontal axis) at two classic EUC
# longitudes — 140°W and 110°W (idealised basin assumed to span 120°E→80°W):
plot_meridional_last100(filepath, output_dir; window_days = 100, depth_max = 400.0,
                        xloc = pacific_lon_to_x(-140.0), label = "≈140°W", tag = "140W")
plot_meridional_last100(filepath, output_dir; window_days = 100, depth_max = 400.0,
                        xloc = pacific_lon_to_x(-110.0), label = "≈110°W", tag = "110W")

# Surface fields (last-100-day mean):
plot_surface_last100(filepath, output_dir; window_days = 100)

# Surface heat budget — reconstructs the restoring flux from the saved SST, so its
# λ_T / target params MUST match toy_model.jl (λ_T doubled to 4e-5 on 2026-10-20;
# the function's own default is still 2e-5).
plot_surface_heat_budget(filepath, output_dir; λ_T = 4e-5)



println("Finished run.")
println("Outputs saved in: ", output_dir)

