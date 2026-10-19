#!/bin/bash
# ============================================================================
# PBS wrapper to run run_toy_model.jl (full run + plots) on Gadi (NCI).
# Submit from the toy_model/ directory:   cd toy_model && qsub run_toy_model.sh
# Runs whatever run_toy_model.jl is set to (stop_time / switches live there),
# CPU / multithreaded on the `normal` queue.
# ============================================================================

#PBS -P e14
#PBS -q normal
#PBS -l ncpus=12
#PBS -l mem=48GB
#PBS -l walltime=10:00:00
#PBS -l storage=scratch/e14+gdata/e14      # repo+outputs on scratch/e14, julia+depot on gdata/e14
#PBS -l wd                                  # start in the dir qsub was run from
#PBS -N toy_model
#PBS -j oe
#PBS -o toy_model.log

set -e

# --- Julia environment (fresh personal install in project space; see SETUP below) ---
# Depot MUST match the one used at setup (keeps packages/artifacts off home's 10 GB):
export JULIA_DEPOT_PATH=/g/data/e14/$USER/.julia    # = /g/data/e14/ic0706/.julia (your existing depot)
export PATH=/g/data/e14/$USER/julia-1.10/bin:$PATH  # <<< adjust if you install julia elsewhere

# --project=. uses the Project.toml/Manifest.toml created in toy_model/ at setup.
julia --project=. -t $PBS_NCPUS run_toy_model.jl

# ============================================================================
# GPU alternative (for a larger grid later — your tilted-BL runs were CPU). To use it:
#   1. change RectilinearGrid(CPU(), ...) -> RectilinearGrid(GPU(), ...) in toy_model.jl
#   2. replace the #PBS resource lines above with:
#        #PBS -q gpuvolta
#        #PBS -l ngpus=1
#        #PBS -l ncpus=12
#        #PBS -l mem=90GB
#      (gpuvolta ratio: 1 GPU : 12 CPUs : 90 GB per V100)
#   3. run with:   julia run_toy_model.jl     (no -t needed on GPU)
# ============================================================================
