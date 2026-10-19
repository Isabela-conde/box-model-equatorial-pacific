# Running the box model on Gadi (NCI)

One-time Julia setup + how to submit a run. Committed so it isn't lost.

Filled in for project **e14**, user **ic0706** (depot `/g/data/e14/ic0706/.julia`).

> **Your `/g/data/e14/ic0706/.julia` depot probably still exists** (`/g/data` isn't
> purged like `/scratch`), so Oceananigans etc. may already be cached — only the julia
> **binary** needs reinstalling (§1), and `Pkg.add` (§2) should be fast. Check with
> `ls /g/data/e14/ic0706/.julia/packages` before assuming a full re-download.

## 0. Key gotchas
- **Login nodes have internet; compute nodes do NOT.** All package downloads
  (`Pkg.add`, `Pkg.instantiate`) must be done on the **login node** before you `qsub`.
- **Depot on /g/data, not home.** Julia packages + artifacts (CairoMakie, CUDA)
  overflow home's 10 GB. Point `JULIA_DEPOT_PATH` at project space.
- Keep the **repo + outputs on /scratch** (runs/ can get big); **Julia + depot on /g/data**.

## 1. Install Julia (login node, one time)
```bash
PROJ=e14
cd /g/data/$PROJ/$USER
# grab the current 1.10.x linux-x86_64 tarball from https://julialang.org/downloads/
wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.9-linux-x86_64.tar.gz
tar zxf julia-1.10.9-linux-x86_64.tar.gz
ln -s julia-1.10.9 julia-1.10          # stable name so run_toy_model.sh's PATH survives patch bumps
```

## 2. Build the project environment (login node, one time)
```bash
PROJ=e14
export PATH=/g/data/$PROJ/$USER/julia-1.10/bin:$PATH
export JULIA_DEPOT_PATH=/g/data/$PROJ/$USER/.julia
mkdir -p $JULIA_DEPOT_PATH

# clone the repo into scratch (outputs go to toy_model/runs/)
cd /scratch/$PROJ/$USER
git clone <repo-url> && cd */toy_model     # cd into the toy_model dir

# download deps now (login node = internet). Skip precompile here (login-node time
# limits) — the first job will precompile on its own.
export JULIA_PKG_PRECOMPILE_AUTO=0
julia --project=. -e 'using Pkg; Pkg.add(["Oceananigans","JLD2","SeawaterPolynomials","CairoMakie"]); Pkg.instantiate()'
```

## 3. Submit the run
```bash
# edit run_toy_model.sh: fill the <<< spots (PROJECT code, and check the julia paths)
cd /scratch/$PROJ/$USER/*/toy_model
qsub run_toy_model.sh
```
- The **first job spends ~10–20 min precompiling** (Oceananigans + CairoMakie) before
  the model starts — normal, only happens once per depot.
- CairoMakie is headless, so the end-of-run plots work on the compute node.
- Monitor: `qstat -u $USER`; log streams to `toy_model.log`.

## Notes
- CPU / `normal` queue is the default here (small grid, 140×40×60). For GPU, see the
  commented block at the bottom of `run_toy_model.sh` (needs `CPU()`→`GPU()` in toy_model.jl
  and the `gpuvolta` queue).
- Current run_toy_model.jl config: relative wind ON, poleward sink ON, 400 days.
