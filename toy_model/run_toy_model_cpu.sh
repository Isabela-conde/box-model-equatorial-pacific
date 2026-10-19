#!/bin/bash
# ============================================================================
# CPU MIRROR of run_toy_model.sh.  Runs run_toy_model_cpu.jl (same model, same
# run + plots, but a CPU() grid) on the `normal` queue.  Use it as a check, or as a
# backup when gpuvolta is slow to schedule.  The GPU setup is untouched.
# Submit from the toy_model/ directory:   cd toy_model && qsub run_toy_model_cpu.sh
#
# NOTE: at the current 600x150x60 grid a CPU run is SLOW (~16x the old 140x40 grid;
# roughly a day+ per 30 model days on 12 cores).  For a quick check, lower Nx/Ny in
# toy_model_cpu.jl (e.g. Nx=140, Ny=60 = the 55 km grid) — the GPU file is unaffected.
# ============================================================================

#PBS -P e14
#PBS -q normal
#PBS -l ncpus=12
#PBS -l mem=48GB
#PBS -l walltime=48:00:00
#PBS -l storage=scratch/e14+gdata/e14
#PBS -l wd
#PBS -N toy_model_cpu
#PBS -j oe
#PBS -o toy_model_cpu.log

set -e

# --- Julia environment (same depot/binary as the GPU run) ---
export JULIA_DEPOT_PATH=/g/data/e14/$USER/.julia
export PATH=/g/data/e14/$USER/julia-1.10/bin:$PATH

# The cuda module IS needed even on CPU: Oceananigans always loads CUDA.jl, and with
# the local-toolkit preference its jlls (CUDA_Compiler_jll etc.) must find the toolkit
# to initialise — otherwise they try to download an artifact, which fails offline.
# No GPU is used; this just lets CUDA.jl load.
module load cuda/12.6.2
export JULIA_CUDA_USE_COMPAT=false

# LIVE LOG: PBS only delivers toy_model_cpu.log when the job ENDS, so julia's output
# is redirected to run_cpu_<jobid>.out on /scratch — `tail -f` it to watch progress.
# -t uses all 12 cores (Oceananigans CPU is multithreaded).
julia --project=. -t $PBS_NCPUS run_toy_model_cpu.jl > run_cpu_${PBS_JOBID%%.*}.out 2>&1
