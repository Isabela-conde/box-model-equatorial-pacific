# Re-plot an EXISTING run without re-running the model.
#
# All plotting functions now live in toy_plots.jl (single source of truth) and
# are also called automatically at the end of run_toy_model.jl.  This script is
# just for re-plotting an old run: point `output_dir` at it and run the file.

include("toy_plots.jl")   # resolves relative to this file's directory

################################################
#### Point at a run directory and plot ####
################################################

repo_root  = dirname(@__DIR__)
output_dir = joinpath(repo_root, "runs", "toy_model_2026-06-05_234636")  # ← change to your run
filepath   = joinpath(output_dir, "toy_model.jld2")

# ȳ-mean and equatorial-band cross sections:
plot_ymean_last100(filepath, output_dir; window_days = 100, depth_max = 400.0)
plot_ymean_last100(filepath, output_dir; window_days = 100, depth_max = 400.0, y_band = 250e3)

# Surface fields (last-100-day mean):
plot_surface_last100(filepath, output_dir; window_days = 100)

# Surface heat budget — PARAMS MUST MATCH THE RUN.
# Defaults match the current model (λ_T=7e-6, evaporation, 70 W/m² solar):
plot_surface_heat_budget(filepath, output_dir)

# Momentum spin-down diagnostic, if the file is present:
momentum_file = joinpath(output_dir, "toy_model_momentum.jld2")
isfile(momentum_file) && plot_momentum_diagnostics(momentum_file, output_dir)

println("Done — plots saved in ", output_dir)
