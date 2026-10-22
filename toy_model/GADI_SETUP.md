# Running the box model on Gadi (NCI) — GPU

One-time setup plus the day-to-day workflow. Filled in for project **e14**, user
**ic0706**: Julia at `/g/data/e14/ic0706/julia-1.10`, depot `/g/data/e14/ic0706/.julia`,
repo at `/scratch/e14/ic0706/box-model-equatorial-pacific`. Runs use one V100 on the
`gpuvolta` queue (Oceananigans 0.96.35, CUDA.jl 5.8, CUDA toolkit 12.6); a model year
takes ~1.3 h of GPU wall time.

## 0. Things that bite
- **Login nodes have internet; compute nodes do not.** Every `Pkg.add` / `Pkg.instantiate`
  / artifact download happens on the login node, before `qsub`.
- **Depot and Julia binary on `/g/data`, repo and outputs on `/scratch`.** The depot
  (CairoMakie, CUDA artifacts) overflows home's 10 GB; `runs/` gets big (a 5-year
  branch is ~30 GB) and belongs on scratch.
- **V100 needs CUDA 12, not 13.** CUDA.jl is pinned to 5.8 and told to use the module's
  toolkit (`LocalPreferences.toml`); without the pin the compiler artifact drops `sm_70`
  and the run dies with `ptxas fatal: Value 'sm_70' is not defined`.
- **`#PBS` lines take no trailing comments.** `qsub` rejects the whole directive.
- **PBS spools stdout until the job ends.** `run_toy_model.sh` redirects Julia's output to
  `run_<jobid>.out` in `toy_model/` so it can be followed live with `tail -f`.
- Commits and pushes happen on the laptop (GitHub Desktop); Gadi only pulls.

## 1. Julia binary (login node, once)
```bash
cd /g/data/e14/ic0706
wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.9-linux-x86_64.tar.gz
tar zxf julia-1.10.9-linux-x86_64.tar.gz
ln -s julia-1.10.9 julia-1.10          # stable name used by run_toy_model.sh
```

## 2. Environment for every interactive session (login node)
```bash
cd /scratch/e14/ic0706/box-model-equatorial-pacific/toy_model
export JULIA_DEPOT_PATH=/g/data/e14/ic0706/.julia
export PATH=/g/data/e14/ic0706/julia-1.10/bin:$PATH
module load cuda/12.6.2
export JULIA_CUDA_USE_COMPAT=false
```
Oceananigans loads CUDA.jl unconditionally, so the cuda module is needed even for
analysis scripts that never touch the GPU.

## 3. Project environment (login node, once)
```bash
# (environment from §2 first)
export JULIA_PKG_PRECOMPILE_AUTO=0        # precompile happens in the first job instead
julia --project=. -e 'using Pkg; Pkg.add(name="CUDA", version="5.8"); Pkg.pin("CUDA");
                      Pkg.add(["Oceananigans","JLD2","SeawaterPolynomials","CairoMakie"]); Pkg.instantiate()'
julia --project=. -e 'using CUDA; CUDA.set_runtime_version!(v"12.6"; local_toolkit=true)'
```
This writes `Project.toml`, `Manifest.toml` and `LocalPreferences.toml` in `toy_model/`
(all untracked — they are Gadi-specific). Check:
```bash
grep -A2 'name = "CUDA"' Manifest.toml      # version = "5.8.x"
cat LocalPreferences.toml                     # [CUDA_Runtime_jll] local = "true", version = "12.6"
julia --project=. -e 'using CUDA; println("CUDA ok"); include("toy_model.jl"); println("scripts ok")'
```

## 4. Pull the latest code
```bash
cd /scratch/e14/ic0706/box-model-equatorial-pacific && git fetch origin && git reset --hard origin/main && git log -1 --oneline
```
`reset --hard` discards local edits to tracked files but leaves the untracked
`Project.toml` / `Manifest.toml` / `LocalPreferences.toml` / `runs/` alone.
Do not add `git clean`.

## 5. Submit
What a job runs is chosen at submit time with the `BRANCH` variable, read by
`run_toy_model.jl` (see its header for the full list); `BRANCH_YEARS` sets the branch
length. Several branches can sit in the queue together.
```bash
cd /scratch/e14/ic0706/box-model-equatorial-pacific/toy_model

# 5-year spin-up from rest (≈3.5 h incl. plots)
qsub run_toy_model.sh

# 1-year branches from the spin-up checkpoint (≈35 min each)
qsub -v BRANCH=control -N control -o control.log -l walltime=02:00:00 run_toy_model.sh
qsub -v BRANCH=warmA2  -N warmA2  -o warmA2.log  -l walltime=02:00:00 run_toy_model.sh

# 5-year branches (≈3 h each)
for b in control warmA1 warmA2 warmA4 warmB2 warmC2; do
  qsub -v BRANCH=$b,BRANCH_YEARS=5 -N ${b}_5y -o ${b}_5y.log -l walltime=06:00:00 run_toy_model.sh
done
```
Branches start from `spinup_ckpt` in `run_toy_model.jl` (the production spin-up's
day-1800 checkpoint); after a new spin-up, update that path.

Monitor:
```bash
qstat -swu ic0706                      # your jobs + scheduler comment
tail -f run_<jobid>.out                # live model log (iteration, model time, Δt, max|w|)
```
The first job after a fresh depot spends ~10–20 min precompiling. Each run writes to
`runs/toy_model_<branch>_<timestamp>/`: 3-D time means, 2-day slice files
(`_surface`, `_140W`, `_110W`, `_ymean`), momentum series, checkpoint, plots.

On runs longer than ~5 years the end-of-run plots can exceed the 90 GB memory
request (they load the whole 3-D file); the job is then killed *after* all data files
are written, so nothing is lost. Use `analysis/` for long runs.

## 6. Analysis on Gadi (login node)
```bash
# (environment from §2 first)
julia --project=. analysis/spinup_diagnostics.jl runs/toy_model_spinup_<ts>
julia --project=. analysis/band_section.jl runs/toy_model_<branch>_<ts> 365 5   # small x–z section file
```
`band_section.jl` reduces a run's 3-D file to a few-MB `band_section_5deg.jld2`; fetch
those (and the slice files) to the laptop rather than the 3-D files.

## 7. Download to the laptop
```bash
# slices + plots, no 3-D file or checkpoint
rsync -avh --progress --exclude 'toy_model.jld2' --exclude 'checkpoint_*' \
  ic0706@gadi.nci.org.au:/scratch/e14/ic0706/box-model-equatorial-pacific/toy_model/runs/toy_model_<branch>_<ts> ~/Downloads/
```
Laptop analysis (`analysis/*.jl`) reads the files raw with JLD2, so the laptop's newer
Oceananigans is fine; it cannot open the Gadi files with `FieldTimeSeries`.
