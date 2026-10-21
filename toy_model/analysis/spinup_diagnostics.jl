# =============================================================================
# spinup_diagnostics.jl — equilibrium / budget diagnostics for ONE run folder
# (a spin-up or a branch), from the files run_model writes:
#     toy_model.jld2           3-D time means (T,S,u,v,w,dTdz)
#     toy_model_surface.jld2   top-cell slices, every 2 days
#     toy_model_140W.jld2, toy_model_110W.jld2   meridional sections, every 2 days
#     toy_model_ymean.jld2     y-mean T, u (x–z), every 2 days
#
# Usage — three ways to say which run:
#   1. edit RUNDIR below (LAPTOP default: a run folder downloaded to ~/Downloads) and run
#          julia analysis/spinup_diagnostics.jl              (or include() it from a REPL)
#   2. give the folder on the command line (overrides RUNDIR):
#          julia --project=. analysis/spinup_diagnostics.jl runs/toy_model_spinup_2026-10-21_140641
#   3. give the folder and the averaging window in days at the END of the run (default 100):
#          julia --project=. analysis/spinup_diagnostics.jl runs/toy_model_warmA2_<ts> 100
#
# Writes into OUTDIR (default: the run folder itself, so a laptop run saves next to the
# downloaded data):
#     diag_timeseries.png     EUC core/depth, surface u, z20, SST, SSS, 2Δx noise vs time
#     diag_heat_budget.png    dHC/dt vs surface restoring + poleward sink (+ residual)
#     diag_euc_structure.png  u, T, N², Ri profiles at 140W/110W over the window
#     diag_eqband_section.png equatorial-band T and u sections over the window
# and prints a summary table (compare with the TAO numbers printed alongside).
#
# The files are read RAW with JLD2 (timeseries/<var>/<iteration>, 4-cell halos), so this
# works with any Oceananigans version — FieldTimeSeries from 0.107 cannot open 0.96 files.
# Grid and forcing parameters are the production values of toy_model.jl (2026-10-21);
# warming branches are recognised from the folder name (warmA1/A2/A4, warmB2, warmC2) so
# the budget uses the shifted restoring targets.
# =============================================================================
using JLD2, Statistics, Printf, CairoMakie
using SeawaterPolynomials, SeawaterPolynomials.TEOS10

# ---------------------------------------------------------------------------
# WHICH RUN (laptop defaults — edit these two lines; command-line arguments override them)
# ---------------------------------------------------------------------------
RUNDIR = expanduser("~/Downloads/toy_model_spinup_2026-10-21_140641")   # downloaded run folder
WINDOW = 100.0                                                           # averaging window (days)
OUTDIR = nothing            # nothing → save the figures in RUNDIR; or e.g. expanduser("~/Desktop/diag")

rundir  = abspath(length(ARGS) ≥ 1 ? ARGS[1] : RUNDIR)
window  = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : WINDOW
outdir  = OUTDIR === nothing ? rundir : abspath(OUTDIR)
mkpath(outdir)
runname = basename(rundir)
isdir(rundir) || error("run folder not found: $rundir  (edit RUNDIR at the top of the script, or pass the folder as the first argument)")
isfile(joinpath(rundir, "toy_model_surface.jld2")) || error("no toy_model_surface.jld2 in $rundir — is the download complete?")
println("=== spinup_diagnostics: $runname  (window = last $(round(Int, window)) days)")
println("    reading  $rundir\n    writing  $outdir")

# ---------------------------------------------------------------------------
# configuration: production toy_model.jl (2026-10-21).  Edit here if the model changes.
# ---------------------------------------------------------------------------
const Nx, Ny, Nz = 600, 150, 60
const Lx, Ly, Lz = 1.5e7, 3.33e6, 2000.0
const refinement, stretching = 3.5, 18
const H = 4                                   # halo cells in the saved arrays
const ρ₀, cₚ, g = 1026.0, 3991.0, 9.81
const day = 86400.0

# surface restoring (toy_model.jl §4.1) and poleward sink (§5.2)
λ_T = 4e-5;  T_warm0 = 29.5;  T_cold0 = 24.0;  x_front = 0.6;  front_width = 0.18
λ_sink = 1/(150day);  T_surf_ref0 = 22.0;  T_deep_ref = 3.5;  z_ref = -150.0;  d_ref = 120.0
z_sink_mid = -300.0;  d_sink = 75.0;  sponge_width_y = 5e5

# warming branches shift the targets — recognise them from the folder name
ΔT_target = 0.0; ΔT_east = 0.0; ΔT_sink = 0.0
if (m = match(r"warmA(\d+)", runname)) !== nothing
    ΔT_target = parse(Float64, m[1])
elseif occursin("warmB2", runname)
    ΔT_target = 2.0; ΔT_sink = 2.0
elseif occursin("warmC2", runname)
    ΔT_target = 2.0; ΔT_east = 2.0
end
T_warm = T_warm0 + ΔT_target;  T_cold = T_cold0 + ΔT_target + ΔT_east;  T_surf_ref = T_surf_ref0 + ΔT_sink
println("branch offsets: ΔT_target = $ΔT_target  ΔT_east = $ΔT_east  ΔT_sink = $ΔT_sink  → targets $T_warm / $T_cold °C, sink ref $T_surf_ref °C")

# ---------------------------------------------------------------------------
# grid (rebuilt from the same formulas as toy_model.jl)
# ---------------------------------------------------------------------------
hk(k) = (k-1)/Nz;  ζ₀(k) = 1 + (hk(k)-1)/refinement
Σ(k)  = (1 - exp(-stretching*hk(k)))/(1 - exp(-stretching))
zF = [Lz*(ζ₀(k)*Σ(k) - 1) for k in 1:Nz+1];  zC = 0.5 .* (zF[1:end-1] .+ zF[2:end]);  dz = diff(zF);  zI = zF[2:end-1]
xC = [(i-0.5)*Lx/Nx for i in 1:Nx];  xF = [(i-1)*Lx/Nx for i in 1:Nx+1];  yC = [-Ly/2 + (j-0.5)*Ly/Ny for j in 1:Ny]
dA = (Lx/Nx)*(Ly/Ny)
jband = findall(abs.(yC) .<= 250e3)        # equatorial band used by the plots
jcore = findall(abs.(yC) .<= 100e3)        # EUC core band
jeq   = argmin(abs.(yC))
nw    = round(Int, 0.1Nx)                  # "west" / "east" = outer 10 % of the basin
i140  = round(Int, 0.625Nx);  i110 = round(Int, 0.81Nx)
eos   = TEOS10EquationOfState(reference_density = ρ₀)

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
# strip halos on every dimension that has them (sliced / averaged dims have size 1)
strip(A) = A[ntuple(d -> size(A, d) > 1 ? (H+1:size(A, d)-H) : (1:1), ndims(A))...]
iters_times(f, var) = begin
    its = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/$var"]))))
    its, [f["timeseries/t/$k"] for k in its] ./ day
end
function z20(Tcol)                          # depth of the 20 °C isotherm (linear interpolation)
    k = findlast(Tcol .< 20.0)
    (k === nothing || k == length(zC)) && return NaN
    zC[k] + (20.0 - Tcol[k]) * (zC[k+1] - zC[k]) / (Tcol[k+1] - Tcol[k])
end
d2rms(v, idx) = sqrt(mean((v[1:end-2] .- 2v[2:end-1] .+ v[3:end])[idx].^2))   # 2Δx checkerboard metric
Ttarget(x) = T_cold + (T_warm - T_cold)*0.5*(1 - tanh((x/Lx - x_front)/front_width))
nsmask(y)  = max(max(0, (sponge_width_y - (y + Ly/2))/sponge_width_y)^2, max(0, (y - (Ly/2 - sponge_width_y))/sponge_width_y)^2)
Tref(z)    = T_deep_ref + (T_surf_ref - T_deep_ref)*0.5*(1 + tanh((z - z_ref)/d_ref))
zmask(z)   = 0.5*(1 + tanh((z - z_sink_mid)/d_sink))
mask_y = nsmask.(yC);  Tref_z = Tref.(zC);  zmask_z = zmask.(zC);  Ttar_x = Ttarget.(xC)

# =============================================================================
# 1. time series from the 2-day slice files
# =============================================================================
println("\n--- 1. time series (slice files)")
sec = Dict{String,Any}()
for (name, path) in (("140W", "toy_model_140W.jld2"), ("110W", "toy_model_110W.jld2"))
    f = jldopen(joinpath(rundir, path), "r"); its, t = iters_times(f, "u")
    euc = Float64[]; zeuc = Float64[]; usurf = Float64[]; ucore_surf = Float64[]
    for k in its
        u = strip(f["timeseries/u/$k"])[1, :, :]
        ub = dropdims(mean(u[jband, :]; dims = 1); dims = 1); uc = dropdims(mean(u[jcore, :]; dims = 1); dims = 1)
        kz = findall(-350 .<= zC .<= -20); kk = kz[argmax(ub[kz])]
        push!(euc, ub[kk]); push!(zeuc, zC[kk]); push!(usurf, ub[end]); push!(ucore_surf, uc[end])
    end
    close(f); sec[name] = (; t, euc, zeuc, usurf, ucore_surf)
end
fy = jldopen(joinpath(rundir, "toy_model_ymean.jld2"), "r"); itsy, ty = iters_times(fy, "T")
z20w = Float64[]; z20c = Float64[]; z20e = Float64[]; z20rms = Float64[]
for k in itsy
    T = strip(fy["timeseries/T/$k"])[:, 1, :]; zz = [z20(T[i, :]) for i in 1:Nx]
    push!(z20w, mean(filter(!isnan, zz[1:nw]))); push!(z20e, mean(filter(!isnan, zz[Nx-nw+1:Nx]))); push!(z20c, zz[i140])
    d = filter(!isnan, zz[1:end-2] .- 2zz[2:end-1] .+ zz[3:end]); push!(z20rms, sqrt(mean(d.^2)))
end
close(fy)
fs = jldopen(joinpath(rundir, "toy_model_surface.jld2"), "r"); itss, ts = iters_times(fs, "T")
sstw = Float64[]; sste = Float64[]; sstmin = Float64[]; sstmax = Float64[]; sss = Float64[]; sssmin = Float64[]; sssmax = Float64[]
sst2dx = Float64[]; u2dx = Float64[]; Qsurf = Float64[]
ie = round(Int, 0.6Nx):Nx-2
for k in itss
    T = strip(fs["timeseries/T/$k"])[:, :, 1]; S = strip(fs["timeseries/S/$k"])[:, :, 1]; u = strip(fs["timeseries/u/$k"])[1:Nx, :, 1]
    push!(sstw, mean(T[1:nw, jband])); push!(sste, mean(T[Nx-nw+1:Nx, jband])); push!(sstmin, minimum(T)); push!(sstmax, maximum(T))
    push!(sss, mean(S)); push!(sssmin, minimum(S)); push!(sssmax, maximum(S))
    push!(sst2dx, d2rms(T[:, jeq], ie)); push!(u2dx, d2rms(u[:, jeq], ie))
    push!(Qsurf, -ρ₀*cₚ*λ_T*dA*sum(T .- Ttar_x))          # restoring heat flux INTO the ocean (W)
end
close(fs)
tend = ts[end]; wsel = ts .>= tend - window
@printf "run covers days %.0f–%.0f; window = days %.0f–%.0f\n" ts[1] tend tend-window tend
for name in ("140W", "110W")
    s = sec[name]; w = s.t .>= tend - window
    @printf "  %s: EUC core (band ±250 km) %.2f ± %.2f m/s at %.0f m | surface u: band %+.2f, core band ±100 km %+.2f\n" name mean(s.euc[w]) std(s.euc[w]) mean(s.zeuc[w]) mean(s.usurf[w]) mean(s.ucore_surf[w])
end
wy = ty .>= tend - window
@printf "  z20 west / 140W / east: %.0f / %.0f / %.0f m   (drift over the window: %+.1f / %+.1f m)\n" mean(z20w[wy]) mean(z20c[wy]) mean(z20e[wy]) z20w[end]-z20w[findfirst(wy)] z20e[end]-z20e[findfirst(wy)]
@printf "  SST west / east band: %.2f / %.2f °C (min %.2f, max %.2f) | SSS mean %.3f (%.3f–%.3f)\n" mean(sstw[wsel]) mean(sste[wsel]) mean(sstmin[wsel]) mean(sstmax[wsel]) mean(sss[wsel]) mean(sssmin[wsel]) mean(sssmax[wsel])
@printf "  2Δx noise along the equator, east 40%%: SST %.3f K, surface u %.3f m/s, y-mean z20 %.2f m\n" mean(sst2dx[wsel]) mean(u2dx[wsel]) mean(z20rms[wy])

fig = Figure(size = (1500, 1100))
Label(fig[0, 1:2], "$runname — time series from the 2-day slice files", fontsize = 18, tellwidth = false)
a = Axis(fig[1, 1], title = "EUC core speed (max band-mean u, 20–350 m)", ylabel = "m/s", xlabel = "day")
lines!(a, sec["140W"].t, sec["140W"].euc, label = "140W"); lines!(a, sec["110W"].t, sec["110W"].euc, label = "110W"); axislegend(a, position = :rb)
a = Axis(fig[1, 2], title = "EUC core depth", ylabel = "m", xlabel = "day")
lines!(a, sec["140W"].t, sec["140W"].zeuc, label = "140W"); lines!(a, sec["110W"].t, sec["110W"].zeuc, label = "110W"); axislegend(a, position = :rb)
a = Axis(fig[2, 1], title = "surface u on the equator (band ±250 km)", ylabel = "m/s", xlabel = "day")
lines!(a, sec["140W"].t, sec["140W"].usurf, label = "140W"); lines!(a, sec["110W"].t, sec["110W"].usurf, label = "110W")
hlines!(a, [0.0], color = :black, linestyle = :dash); axislegend(a, position = :rb)
a = Axis(fig[2, 2], title = "z(20 °C), y-mean: west 10 % / 140W / east 10 %", ylabel = "m", xlabel = "day")
lines!(a, ty, z20w, label = "west"); lines!(a, ty, z20c, label = "140W"); lines!(a, ty, z20e, label = "east"); axislegend(a, position = :rb)
a = Axis(fig[3, 1], title = "SST: band west / east, basin min / max  |  SSS mean (right axis)", ylabel = "°C", xlabel = "day")
lines!(a, ts, sstw, label = "SST west"); lines!(a, ts, sste, label = "SST east")
lines!(a, ts, sstmin, color = :gray40, linestyle = :dash, label = "min"); lines!(a, ts, sstmax, color = :gray40, linestyle = :dot, label = "max")
axislegend(a, position = :rt)
a2 = Axis(fig[3, 1], yaxisposition = :right, ylabel = "SSS", yticklabelcolor = :purple); hidespines!(a2); hidexdecorations!(a2)
lines!(a2, ts, sss, color = :purple)
a = Axis(fig[3, 2], title = "2Δx checkerboard metric (rms 2nd-diff along the equator, east 40 %)", xlabel = "day", yscale = log10)
lines!(a, ts, max.(sst2dx, 1e-5), label = "SST (K)"); lines!(a, ts, max.(u2dx, 1e-5), label = "surface u (m/s)"); lines!(a, ty, max.(z20rms, 1e-5), label = "y-mean z20 (m)")
axislegend(a, position = :rt)
save(joinpath(rundir, "diag_timeseries.png"), fig); println("  wrote diag_timeseries.png")

# =============================================================================
# 2. heat budget: dHC/dt (3-D means) vs surface restoring (slices) + poleward sink (3-D means)
# =============================================================================
println("\n--- 2. heat budget (3-D time means + surface slices)")
f3 = jldopen(joinpath(rundir, "toy_model.jld2"), "r"); its3, t3 = iters_times(f3, "T")
HC = Float64[]; SINK = Float64[]; Tlast = nothing; Slast = nothing; nwin = 0
for (n, k) in enumerate(its3)
    T = strip(f3["timeseries/T/$k"])
    push!(HC, ρ₀*cₚ*dA*sum(T .* reshape(dz, 1, 1, Nz)))
    s = -λ_sink .* reshape(mask_y, 1, Ny, 1) .* reshape(zmask_z, 1, 1, Nz) .* (T .- reshape(Tref_z, 1, 1, Nz))
    push!(SINK, ρ₀*cₚ*dA*sum(s .* reshape(dz, 1, 1, Nz)))
    if t3[k == its3[end] ? n : n] >= tend - window              # accumulate the window mean of the 3-D fields
        global Tlast = Tlast === nothing ? T : Tlast .+ T
        global Slast = Slast === nothing ? strip(f3["timeseries/S/$k"]) : Slast .+ strip(f3["timeseries/S/$k"])
        global nwin += 1
    end
end
close(f3); Tlast ./= nwin; Slast ./= nwin
println("  3-D file: $(length(its3)) snapshots, days $(t3[1])–$(t3[end]); window uses the last $nwin")
dHC = Float64[]; qs = Float64[]; sk = Float64[]; tb = Float64[]
for n in 2:length(t3)
    a, b = t3[n-1], t3[n]
    push!(tb, 0.5(a+b)); push!(dHC, (HC[n]-HC[n-1])/((b-a)*day))
    push!(qs, mean(Qsurf[(ts .> a) .& (ts .<= b)])); push!(sk, 0.5(SINK[n]+SINK[n-1]))
end
res = dHC .- qs .- sk
wb = tb .>= tend - window
@printf "  window means: dHC/dt %+.2f PW | surface restoring %+.2f PW | sink %+.2f PW | residual %+.2f PW (= %+.1f W/m²)\n" mean(dHC[wb])/1e15 mean(qs[wb])/1e15 mean(sk[wb])/1e15 mean(res[wb])/1e15 mean(res[wb])/(Lx*Ly)
println("  (residual ≈ 0 ⇒ the model conserves heat; dHC/dt ≈ 0 ⇒ equilibrated)")
fig = Figure(size = (1200, 700))
Label(fig[0, 1], "$runname — basin heat budget (PW, + into the ocean)", fontsize = 18, tellwidth = false)
a = Axis(fig[1, 1], xlabel = "day", ylabel = "PW")
lines!(a, tb, dHC ./ 1e15, color = :purple, linewidth = 3, label = "dHC/dt")
lines!(a, tb, qs ./ 1e15, color = :firebrick, label = "surface restoring")
lines!(a, tb, sk ./ 1e15, color = :royalblue, label = "poleward sink")
lines!(a, tb, res ./ 1e15, color = :gray40, linestyle = :dash, label = "residual (non-conservation)")
hlines!(a, [0.0], color = :black, linewidth = 0.8); axislegend(a, position = :rt)
save(joinpath(rundir, "diag_heat_budget.png"), fig); println("  wrote diag_heat_budget.png")

# =============================================================================
# 3. EUC vertical structure and Richardson number at 140W / 110W (window mean, |y| ≤ 100 km)
# =============================================================================
println("\n--- 3. EUC structure over the window (|y| ≤ 100 km)")
prof = Dict{String,Any}()
for (name, path) in (("140W", "toy_model_140W.jld2"), ("110W", "toy_model_110W.jld2"))
    f = jldopen(joinpath(rundir, path), "r"); its, t = iters_times(f, "u"); sel = its[t .>= tend - window]
    p(v) = mean([dropdims(mean(strip(f["timeseries/$v/$k"])[1, jcore, :]; dims = 1); dims = 1) for k in sel])
    u = p("u"); v = p("v"); T = p("T"); S = p("S"); close(f)
    ρ = [SeawaterPolynomials.ρ′(T[k], S[k], zC[k], eos) for k in 1:Nz]
    N2 = -(g/ρ₀) .* diff(ρ) ./ diff(zC); Sh2 = (diff(u) ./ diff(zC)).^2 .+ (diff(v) ./ diff(zC)).^2; Ri = N2 ./ max.(Sh2, 1e-12)
    prof[name] = (; u, T, N2, Ri)
    kc = argmax(u); k0 = findlast(u .< 0)
    @printf "  %s: surface u %+.2f | u at 25/50 m %+.2f/%+.2f | core %+.2f at %.0f m | Ri at 30 m %.2f | N² at 10/50 m %.1e/%.1e\n" name u[end] u[argmin(abs.(zC .+ 25))] u[argmin(abs.(zC .+ 50))] u[kc] zC[kc] Ri[argmin(abs.(zI .+ 30))] N2[argmin(abs.(zI .+ 10))] N2[argmin(abs.(zI .+ 50))]
end
println("  TAO/Johnson at 140W for reference: surface −0.2…−0.5, zero crossing 40–60 m, core +1.0–1.2 at 110–130 m, Ri 0.25–0.5 in the shear layer")
fig = Figure(size = (1500, 500))
Label(fig[0, 1:4], "$runname — equatorial profiles (|y| ≤ 100 km), last $(round(Int, window)) days", fontsize = 18, tellwidth = false)
for (ci, name) in enumerate(("140W", "110W")); p = prof[name]
    a = Axis(fig[1, 2ci-1], title = "$name  u (m/s)", ylabel = "z (m)"); lines!(a, p.u, zC); vlines!(a, [0.0], color = :black, linestyle = :dash); ylims!(a, -300, 0)
    a = Axis(fig[1, 2ci], title = "$name  Ri (log) and T"); lines!(a, max.(p.Ri, 1e-2), zI, label = "Ri"); vlines!(a, [0.25], color = :black, linestyle = :dash)
    a.xscale = log10; ylims!(a, -300, 0)
    a2 = Axis(fig[1, 2ci], xaxisposition = :top, xlabel = "T (°C)", xticklabelcolor = :firebrick); hidespines!(a2); hideydecorations!(a2)
    lines!(a2, p.T, zC, color = :firebrick); ylims!(a2, -300, 0)
end
save(joinpath(rundir, "diag_euc_structure.png"), fig); println("  wrote diag_euc_structure.png")

# =============================================================================
# 4. equatorial-band section (window mean of the 3-D means) + zonal pressure-gradient proxy
# =============================================================================
println("\n--- 4. equatorial-band section and pressure-gradient proxy (window mean of the 3-D file)")
Teq = dropdims(mean(Tlast[:, jband, :]; dims = 2); dims = 2)
f3 = jldopen(joinpath(rundir, "toy_model.jld2"), "r")
sel3 = its3[t3 .>= tend - window]
ueq = mean([dropdims(mean(strip(f3["timeseries/u/$k"])[:, jband, :]; dims = 2); dims = 2) for k in sel3]); close(f3)
ρa = [mean(SeawaterPolynomials.ρ′(Tlast[i, j, k], Slast[i, j, k], zC[k], eos) for j in jcore) for i in 1:Nx, k in 1:Nz]
ηs = [-(1/ρ₀)*sum(ρa[i, :] .* dz) for i in 1:Nx]                       # steric height rel. 2000 m
iw = argmin(abs.(xC .- 2e6)); ie2 = argmin(abs.(xC .- 13e6))
pgf = ρ₀*g*(ηs[iw] - ηs[ie2])/(xC[ie2] - xC[iw])
@printf "  steric height W−E (2000→13000 km, equator) = %.2f m (real Pacific ≈ 0.55) → surface PGF %.2e N/m³\n" ηs[iw]-ηs[ie2] pgf
kz350 = findall(zC .>= -350); kmax = kz350[argmax(maximum(ueq[:, kz350]; dims = 1)[:])]
@printf "  band-mean EUC max over the section: %+.2f m/s at %.0f m; z20 along the equator W/E: %.0f / %.0f m\n" maximum(ueq[:, kz350]) zC[kmax] z20(Teq[nw, :]) z20(Teq[Nx-nw, :])
fig = Figure(size = (1400, 850))
Label(fig[0, 1:2], "$runname — equatorial band (|y| ≤ 250 km), last $(round(Int, window)) days", fontsize = 18, tellwidth = false)
a = Axis(fig[1, 1], title = "T (°C), white = 20 °C", ylabel = "z (m)"); hm = heatmap!(a, xC ./ 1e3, zC, Teq, colormap = :thermal, colorrange = (10, 30))
contour!(a, xC ./ 1e3, zC, Teq, levels = [20.0], color = :white, linewidth = 2); ylims!(a, -400, 0); Colorbar(fig[1, 2], hm)
xu = size(ueq, 1) == Nx + 1 ? xF : xC
a = Axis(fig[2, 1], title = "u (m/s), contours every 0.1", xlabel = "x (km, W→E)", ylabel = "z (m)")
hm = heatmap!(a, xu ./ 1e3, zC, ueq, colormap = :balance, colorrange = (-0.8, 0.8))
contour!(a, xu ./ 1e3, zC, ueq, levels = -0.7:0.1:0.7, color = :black, linewidth = 0.5); ylims!(a, -400, 0); Colorbar(fig[2, 2], hm)
save(joinpath(rundir, "diag_eqband_section.png"), fig); println("  wrote diag_eqband_section.png")
println("\n=== done: $runname")
