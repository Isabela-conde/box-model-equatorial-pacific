#!/bin/bash
# ============================================================================
# PBS wrapper to run run_toy_model.jl (full run + plots) on Gadi (NCI).
# Submit from the toy_model/ directory:   cd toy_model && qsub run_toy_model.sh
# Runs whatever run_toy_model.jl is set to (stop_time / switches live there),
# Single V100 GPU on the `gpuvolta` queue (per the Oceananigans-on-Gadi guide, #4584).
# ============================================================================

#PBS -P e14
#PBS -q gpuvolta
#PBS -l ngpus=1
#PBS -l ncpus=12
#PBS -l mem=90GB
#PBS -l walltime=04:00:00
# storage: repo+outputs on scratch/e14, julia+depot on gdata/e14
#PBS -l storage=scratch/e14+gdata/e14
# wd: start in the dir qsub was run from
#PBS -l wd
#PBS -N toy_model
#PBS -j oe
#PBS -o toy_model.log

set -e

# --- Julia environment (fresh personal install in project space; see SETUP below) ---
# Depot MUST match the one used at setup (keeps packages/artifacts off home's 10 GB):
export JULIA_DEPOT_PATH=/g/data/e14/$USER/.julia    # = /g/data/e14/ic0706/.julia (your existing depot)
export PATH=/g/data/e14/$USER/julia-1.10/bin:$PATH  # <<< adjust if you install julia elsewhere

# --- CUDA for the V100 (gpuvolta), from the Gadi guide (#4584):
#   pin CUDA 12.x — V100 (compute cap 7.0) is NOT supported on CUDA 13+ — and
#   JULIA_CUDA_USE_COMPAT=false is the key fix several users needed for V100.
module load cuda/12.6.2
export JULIA_CUDA_USE_COMPAT=false
export JULIA_CUDA_MEMORY_POOL=none

# --project=. uses the Project.toml/Manifest.toml created in toy_model/ at setup.
# LIVE LOG: PBS only delivers toy_model.log when the job ENDS (stdout is spooled on
# the compute node until then). Redirecting here writes run_<jobid>.out to /scratch
# as the job runs, so you can `tail -f` it to watch progress / catch errors early.
# (No -t: threads are irrelevant on the GPU.)
julia --project=. run_toy_model.jl > run_${PBS_JOBID%%.*}.out 2>&1

# ============================================================================
# CPU fallback (the pre-rebuild setup). To go back:
#   1. RectilinearGrid(GPU(), ...) -> RectilinearGrid(CPU(), ...) in toy_model.jl
#   2. #PBS -q normal / -l ncpus=12 / -l mem=48GB  (drop ngpus, the cuda module,
#      and the JULIA_CUDA_* exports)
#   3. run with:   julia --project=. -t $PBS_NCPUS run_toy_model.jl
# ============================================================================
