# =============================================================================
# band_section.jl — reduce a run's 3-D output to a latitude-band-mean (x–z) section
# of T, N² and u over the last `window` days, saved as a SMALL file next to the run.
#
#     julia --project=. analysis/band_section.jl <run_folder> [window_days=365] [half_width_deg=5]
#
# Writes <run_folder>/band_section_<half_width>deg.jld2  (a few MB: x, z, T, N2, u, meta).
# Run this where the 3-D file lives (Gadi login node is fine: one snapshot in memory at a
# time), then download the small file and plot with analysis/plot_band_sections.jl.
# N² is computed per column (TEOS10 density) and THEN band-averaged.
# =============================================================================
using JLD2, Statistics, Printf
using SeawaterPolynomials, SeawaterPolynomials.TEOS10

length(ARGS) ≥ 1 || error("usage: julia band_section.jl <run_folder> [window_days] [half_width_deg]")
rundir = abspath(ARGS[1]); window = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 365.0
halfdeg = length(ARGS) ≥ 3 ? parse(Float64, ARGS[3]) : 5.0
f3path = joinpath(rundir, "toy_model.jld2"); isfile(f3path) || error("no toy_model.jld2 in $rundir")

const Nx, Ny, Nz = 600, 150, 60;  const Lx, Ly, Lz = 1.5e7, 3.33e6, 2000.0
const refinement, stretching = 3.5, 18;  const H = 4;  const g, ρ₀, day = 9.81, 1026.0, 86400.0
hk(k) = (k-1)/Nz;  ζ₀(k) = 1 + (hk(k)-1)/refinement;  Σ(k) = (1 - exp(-stretching*hk(k)))/(1 - exp(-stretching))
zF = [Lz*(ζ₀(k)*Σ(k) - 1) for k in 1:Nz+1];  zC = 0.5 .* (zF[1:end-1] .+ zF[2:end]);  zI = zF[2:end-1]
xC = [(i-0.5)*Lx/Nx for i in 1:Nx];  yC = [-Ly/2 + (j-0.5)*Ly/Ny for j in 1:Ny]
jb = findall(abs.(yC) .<= halfdeg*111e3)
strip(A) = A[ntuple(d -> size(A, d) > 1 ? (H+1:size(A, d)-H) : (1:1), ndims(A))...]
eos = TEOS10EquationOfState(reference_density = ρ₀)

f = jldopen(f3path, "r")
its = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/T"])))); t = [f["timeseries/t/$k"] for k in its] ./ day
sel = its[t .>= t[end] - window]
@printf "%s: %d snapshots in the window (days %.0f–%.0f), band |y| ≤ %.0f° (%d rows)\n" basename(rundir) length(sel) t[end]-window t[end] halfdeg length(jb)
Tm = zeros(Nx, Nz); N2m = zeros(Nx, Nz-1); um = zeros(Nx, Nz)
for (n, k) in enumerate(sel)
    T = strip(f["timeseries/T/$k"])[:, jb, :]; S = strip(f["timeseries/S/$k"])[:, jb, :]; u = strip(f["timeseries/u/$k"])[:, jb, :]
    u = size(u, 1) == Nx + 1 ? 0.5 .* (u[1:end-1, :, :] .+ u[2:end, :, :]) : u          # x-faces → centres
    ρ = [SeawaterPolynomials.ρ′(T[i, j, kk], S[i, j, kk], zC[kk], eos) for i in 1:Nx, j in 1:length(jb), kk in 1:Nz]
    N2 = -(g/ρ₀) .* (ρ[:, :, 2:end] .- ρ[:, :, 1:end-1]) ./ reshape(diff(zC), 1, 1, :)
    Tm .+= dropdims(mean(T; dims = 2); dims = 2); um .+= dropdims(mean(u; dims = 2); dims = 2); N2m .+= dropdims(mean(N2; dims = 2); dims = 2)
    n % 10 == 0 && println("  $n / $(length(sel))")
end
close(f)
Tm ./= length(sel); um ./= length(sel); N2m ./= length(sel)
out = joinpath(rundir, @sprintf("band_section_%.0fdeg.jld2", halfdeg))
jldsave(out; x = xC, zC = zC, zI = zI, T = Tm, N2 = N2m, u = um, window = window, halfdeg = halfdeg,
        tstart = t[end] - window, tend = t[end], run = basename(rundir))
println("wrote $out")
