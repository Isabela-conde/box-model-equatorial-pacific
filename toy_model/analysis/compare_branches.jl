# =============================================================================
# compare_branches.jl — warming branches vs control: EUC / thermocline response and
# the Pedlosky scaling test (H ∝ γ₂^-1/2, U ∝ γ₂^1/4).
#
# Usage — edit the BRANCHES dict below (laptop: downloaded run folders), then
#     julia analysis/compare_branches.jl            (or include() from a REPL)
# Optional first argument = averaging window in days at the END of each branch (default 100).
# Writes into OUTDIR:
#     cmp_profiles.png     u(z), T(z), N²(z) at 140W/110W, all branches, + difference from control
#     cmp_sections.png     equatorial-band Δu and ΔT sections (branch − control)
#     cmp_scaling.png      z20, core depth, core speed, core T vs ΔT and vs γ₂ (log–log, with fitted exponents)
#     cmp_summary.txt      the table printed below
#
# Files are read RAW with JLD2 (any Oceananigans version).  Grid/forcing constants = production
# toy_model.jl (2026-10-21).
# =============================================================================
using JLD2, Statistics, Printf, CairoMakie
using SeawaterPolynomials, SeawaterPolynomials.TEOS10

# ---------------------------------------------------------------------------
# WHICH RUNS (edit here)
# ---------------------------------------------------------------------------
DL = expanduser("~/Downloads")
BRANCHES = [   # (label, folder, ΔT_target, ΔT_east, ΔT_sink)  — control FIRST
    ("control", joinpath(DL, "toy_model_control_2026-10-21_180403"), 0.0, 0.0, 0.0),
    ("A1 +1 K", joinpath(DL, "toy_model_warmA1_2026-10-21_182550"),  1.0, 0.0, 0.0),
    ("A2 +2 K", joinpath(DL, "toy_model_warmA2_2026-10-21_182936"),  2.0, 0.0, 0.0),
    ("A4 +4 K", joinpath(DL, "toy_model_warmA4_2026-10-21_184307"),  4.0, 0.0, 0.0),
    ("B2 +2 K ventilated", joinpath(DL, "toy_model_warmB2_2026-10-21_185704"), 2.0, 0.0, 2.0),
    ("C2 +2/+4 K pattern", joinpath(DL, "toy_model_warmC2_2026-10-21_185857"), 2.0, 2.0, 0.0),
]
WINDOW = length(ARGS) ≥ 1 ? parse(Float64, ARGS[1]) : 100.0
OUTDIR = joinpath(DL, "branch_comparison_2026-10-21")
mkpath(OUTDIR)
cols = [:black, :goldenrod, :firebrick, :purple, :royalblue, :seagreen]
styles = [:solid, :solid, :solid, :solid, :dash, :dot]

# ---------------------------------------------------------------------------
# configuration (production toy_model.jl)
# ---------------------------------------------------------------------------
const Nx, Ny, Nz = 600, 150, 60
const Lx, Ly, Lz = 1.5e7, 3.33e6, 2000.0
const refinement, stretching = 3.5, 18
const H = 4; const ρ₀, cₚ, g = 1026.0, 3991.0, 9.81; const day = 86400.0
const τ_eq = 1.2 * 1.2e-3 * 7.0^2                 # ρₐ cᴰ u₁₀² = 0.071 N/m² (wind is fixed across branches)
hk(k) = (k-1)/Nz;  ζ₀(k) = 1 + (hk(k)-1)/refinement;  Σ(k) = (1 - exp(-stretching*hk(k)))/(1 - exp(-stretching))
zF = [Lz*(ζ₀(k)*Σ(k) - 1) for k in 1:Nz+1];  zC = 0.5 .* (zF[1:end-1] .+ zF[2:end]);  dz = diff(zF);  zI = zF[2:end-1]
xC = [(i-0.5)*Lx/Nx for i in 1:Nx];  xF = [(i-1)*Lx/Nx for i in 1:Nx+1];  yC = [-Ly/2 + (j-0.5)*Ly/Ny for j in 1:Ny]
jband = findall(abs.(yC) .<= 250e3);  jcore = findall(abs.(yC) .<= 100e3);  nw = round(Int, 0.1Nx)
eos = TEOS10EquationOfState(reference_density = ρ₀)
strip(A) = A[ntuple(d -> size(A, d) > 1 ? (H+1:size(A, d)-H) : (1:1), ndims(A))...]
iters_times(f, var) = (its = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/$var"])))); (its, [f["timeseries/t/$k"] for k in its] ./ day))
function zT(Tcol, T0)                       # depth of the T0 isotherm
    k = findlast(Tcol .< T0); (k === nothing || k == length(zC)) && return NaN
    zC[k] + (T0 - Tcol[k])*(zC[k+1]-zC[k])/(Tcol[k+1]-Tcol[k])
end

# ---------------------------------------------------------------------------
# load every branch: window-mean profiles at 140W/110W, band sections, scalar metrics
# ---------------------------------------------------------------------------
R = []
for (lab, dir, dT, dTe, dTs) in BRANCHES
    isdir(dir) || error("missing run folder: $dir")
    r = Dict{Symbol,Any}(:label => lab, :dT => dT, :dTe => dTe, :dTs => dTs)
    for sec in ("140W", "110W")
        f = jldopen(joinpath(dir, "toy_model_$sec.jld2"), "r"); its, t = iters_times(f, "u"); tend = t[end]; sel = its[t .>= tend - WINDOW]
        p(v, jj) = mean([dropdims(mean(strip(f["timeseries/$v/$k"])[1, jj, :]; dims = 1); dims = 1) for k in sel])
        u = p("u", jcore); v = p("v", jcore); T = p("T", jcore); S = p("S", jcore); ub = p("u", jband)
        # core time series (for the ± spread)
        euc_t = [maximum(dropdims(mean(strip(f["timeseries/u/$k"])[1, jband, :]; dims = 1); dims = 1)[zC .>= -350]) for k in sel]
        close(f)
        ρ = [SeawaterPolynomials.ρ′(T[k], S[k], zC[k], eos) for k in 1:Nz]
        N2 = -(g/ρ₀) .* diff(ρ) ./ diff(zC)
        kz = findall(-350 .<= zC .<= -20); kc = kz[argmax(u[kz])]; kb = kz[argmax(ub[kz])]
        r[Symbol(sec)] = (; u, T, N2, ub, core = u[kc], zcore = zC[kc], Tcore = T[kc], core_band = ub[kb], zcore_band = zC[kb],
                           euc_mean = mean(euc_t), euc_std = std(euc_t), usurf = u[end], u25 = u[argmin(abs.(zC .+ 25))],
                           z20 = zT(T, 20.0), zN2max = zI[argmax(N2)], tend)
    end
    # 3-D window mean: band sections, thermocline contrast → γ₂, west/east z20, SST
    f3 = jldopen(joinpath(dir, "toy_model.jld2"), "r"); its3, t3 = iters_times(f3, "T"); tend = t3[end]; sel3 = its3[t3 .>= tend - WINDOW]
    T3 = mean([strip(f3["timeseries/T/$k"]) for k in sel3]); u3 = mean([strip(f3["timeseries/u/$k"]) for k in sel3]); close(f3)
    Teq = dropdims(mean(T3[:, jband, :]; dims = 2); dims = 2); ueq = dropdims(mean(u3[:, jband, :]; dims = 2); dims = 2)
    # reduced gravity across the thermocline: band-mean T at 20 m minus T at 300 m, mid-basin (x = 0.4–0.6 Lx)
    im = round(Int, 0.4Nx):round(Int, 0.6Nx)
    Tup = mean(Teq[im, argmin(abs.(zC .+ 20))]); Tdn = mean(Teq[im, argmin(abs.(zC .+ 300))])
    α = 2.5e-4                                   # thermal expansion at ~25 °C (K⁻¹); salinity contribution small here
    r[:γ₂] = g*α*(Tup - Tdn); r[:ΔTtherm] = Tup - Tdn
    z20x = [zT(Teq[i, :], 20.0) for i in 1:Nx]
    r[:z20w] = mean(filter(!isnan, z20x[1:nw])); r[:z20e] = mean(filter(!isnan, z20x[Nx-nw+1:Nx])); r[:z20m] = mean(filter(!isnan, z20x[im]))
    r[:sstw] = mean(T3[1:nw, jband, end]); r[:sste] = mean(T3[Nx-nw+1:Nx, jband, end])
    kz = findall(zC .>= -350); r[:umax3] = maximum(ueq[:, kz]); r[:Teq] = Teq; r[:ueq] = ueq
    push!(R, r); println("loaded $lab")
end
C = R[1]

# ---------------------------------------------------------------------------
# summary table
# ---------------------------------------------------------------------------
lines = String[]
push!(lines, @sprintf("Window: last %.0f days of each branch (branch day %.0f–%.0f).  Wind fixed, τ_eq = %.3f N/m².", WINDOW, C[Symbol("140W")].tend - WINDOW, C[Symbol("140W")].tend, τ_eq))
push!(lines, @sprintf("%-20s %6s %6s %7s %7s %7s %7s %7s %7s %7s %7s %7s", "branch", "ΔT", "γ₂", "z20mid", "z20W", "z20E", "core140", "zc140", "Tc140", "core110", "zc110", "sstE"))
for r in R
    s = r[Symbol("140W")]; e = r[Symbol("110W")]
    push!(lines, @sprintf("%-20s %+6.1f %6.4f %7.0f %7.0f %7.0f %7.2f %7.0f %7.1f %7.2f %7.0f %7.2f", r[:label], r[:dT], r[:γ₂], r[:z20m], r[:z20w], r[:z20e], s.core, s.zcore, s.Tcore, e.core, e.zcore, r[:sste]))
end
push!(lines, "")
push!(lines, "Response relative to control (branch − control):")
push!(lines, @sprintf("%-20s %8s %8s %8s %8s %8s %8s %8s", "branch", "Δγ₂ %", "Δz20mid", "Δzc140", "Δcore140", "ΔTc140", "Δusurf", "ΔsstE"))
for r in R[2:end]
    s = r[Symbol("140W")]; c = C[Symbol("140W")]
    push!(lines, @sprintf("%-20s %+8.1f %+8.1f %+8.1f %+8.3f %+8.2f %+8.3f %+8.2f", r[:label], 100(r[:γ₂]/C[:γ₂]-1), r[:z20m]-C[:z20m], s.zcore-c.zcore, s.core-c.core, s.Tcore-c.Tcore, s.usurf-c.usurf, r[:sste]-C[:sste]))
end
# ---- scaling fit on the A series (control, A1, A2, A4): log(H), log(U) vs log(γ₂) ----
A = [r for r in R if r[:dTe] == 0 && r[:dTs] == 0]
lg = log.([r[:γ₂] for r in A])
fit(y) = (X = hcat(ones(length(lg)), lg); b = X \ log.(y); b[2])
expH  = fit([-r[:z20m] for r in A]);  expHc = fit([-r[Symbol("140W")].zcore for r in A])
expU  = fit([r[Symbol("140W")].core for r in A]); expUb = fit([r[Symbol("140W")].euc_mean for r in A])
push!(lines, "")
push!(lines, "Pedlosky (1987) scaling test on the A series (wind fixed, so only γ₂ varies):")
push!(lines, @sprintf("  fitted exponent of z20 (mid-basin) vs γ₂:  %+.2f   (inertial prediction −0.50)", expH))
push!(lines, @sprintf("  fitted exponent of EUC core depth vs γ₂:   %+.2f   (inertial prediction −0.50)", expHc))
push!(lines, @sprintf("  fitted exponent of EUC core speed vs γ₂:   %+.2f   (inertial prediction +0.25)", expU))
push!(lines, @sprintf("  fitted exponent of band-mean EUC vs γ₂:    %+.2f   (inertial prediction +0.25)", expUb))
push!(lines, "  (1-year branches: transient response; the thermocline adjusts over years — repeat with BRANCH_YEARS=5)")
foreach(println, lines); open(joinpath(OUTDIR, "cmp_summary.txt"), "w") do io; foreach(l -> println(io, l), lines); end

# ---------------------------------------------------------------------------
# figure 1: profiles
# ---------------------------------------------------------------------------
fig = Figure(size = (1700, 1000))
Label(fig[0, 1:6], "Warming branches vs control — last $(round(Int, WINDOW)) days, |y| ≤ 100 km", fontsize = 18, tellwidth = false)
for (row, sec) in enumerate(("140W", "110W"))
    a1 = Axis(fig[row, 1], title = "$sec  u (m/s)", ylabel = "z (m)"); a2 = Axis(fig[row, 2], title = "$sec  T (°C)")
    a3 = Axis(fig[row, 3], title = "$sec  N² (s⁻²)"); a4 = Axis(fig[row, 4], title = "$sec  Δu vs control (m/s)")
    a5 = Axis(fig[row, 5], title = "$sec  ΔT vs control (K)"); a6 = Axis(fig[row, 6], title = "$sec  ΔN² vs control (s⁻²)")
    for (i, r) in enumerate(R)
        s = r[Symbol(sec)]; c = C[Symbol(sec)]
        lines!(a1, s.u, zC, color = cols[i], linestyle = styles[i], label = r[:label]); lines!(a2, s.T, zC, color = cols[i], linestyle = styles[i]); lines!(a3, s.N2, zI, color = cols[i], linestyle = styles[i])
        i > 1 && (lines!(a4, s.u .- c.u, zC, color = cols[i], linestyle = styles[i]); lines!(a5, s.T .- c.T, zC, color = cols[i], linestyle = styles[i]); lines!(a6, s.N2 .- c.N2, zI, color = cols[i], linestyle = styles[i]))
    end
    vlines!(a1, [0.0], color = :gray, linestyle = :dash); vlines!(a4, [0.0], color = :gray); vlines!(a5, [0.0], color = :gray); vlines!(a6, [0.0], color = :gray)
    for a in (a1, a2, a3, a4, a5, a6); ylims!(a, -300, 0); end
    row == 1 && axislegend(a1, position = :rb, labelsize = 10)
end
save(joinpath(OUTDIR, "cmp_profiles.png"), fig); println("wrote cmp_profiles.png")

# ---------------------------------------------------------------------------
# figure 2: band sections of the response
# ---------------------------------------------------------------------------
fig = Figure(size = (1600, 1500))
Label(fig[0, 1:4], "Equatorial band (|y| ≤ 250 km): branch − control, last $(round(Int, WINDOW)) days.  Black = control 20 °C isotherm, colour = branch", fontsize = 16, tellwidth = false)
xu = size(C[:ueq], 1) == Nx + 1 ? xF : xC
for (i, r) in enumerate(R[2:end])
    row = i
    a = Axis(fig[row, 1], title = "$(r[:label]):  Δu (m/s)", ylabel = "z (m)"); hm = heatmap!(a, xu ./ 1e3, zC, r[:ueq] .- C[:ueq], colormap = :balance, colorrange = (-0.15, 0.15))
    contour!(a, xC ./ 1e3, zC, C[:Teq], levels = [20.0], color = :black, linewidth = 1.5); contour!(a, xC ./ 1e3, zC, r[:Teq], levels = [20.0], color = cols[i+1], linewidth = 1.5); ylims!(a, -400, 0)
    i == 1 && Colorbar(fig[1, 2], hm, label = "Δu (m/s)")
    a = Axis(fig[row, 3], title = "$(r[:label]):  ΔT (K)"); hm2 = heatmap!(a, xC ./ 1e3, zC, r[:Teq] .- C[:Teq], colormap = :balance, colorrange = (-3, 3))
    contour!(a, xC ./ 1e3, zC, C[:Teq], levels = [20.0], color = :black, linewidth = 1.5); contour!(a, xC ./ 1e3, zC, r[:Teq], levels = [20.0], color = cols[i+1], linewidth = 1.5); ylims!(a, -400, 0)
    i == 1 && Colorbar(fig[1, 4], hm2, label = "ΔT (K)")
end
save(joinpath(OUTDIR, "cmp_sections.png"), fig); println("wrote cmp_sections.png")

# ---------------------------------------------------------------------------
# figure 3: scaling
# ---------------------------------------------------------------------------
fig = Figure(size = (1500, 900))
Label(fig[0, 1:3], "Response vs imposed warming (left) and vs reduced gravity γ₂ (right, log–log with the A-series fit)", fontsize = 16, tellwidth = false)
dTs = [r[:dT] + 0.5r[:dTe] for r in R]      # C2 plotted at its basin-mean +3
mk = [i == 1 ? :circle : (r[:dTs] > 0 ? :utriangle : (r[:dTe] > 0 ? :diamond : :circle)) for (i, r) in enumerate(R)]
function panel(pos, title, ylab, y; logx = false, pred = nothing, fitexp = nothing)
    a = Axis(fig[pos...], title = title, ylabel = ylab, xlabel = logx ? "γ₂ (m s⁻²)" : "imposed ΔT (K)", xscale = logx ? log10 : identity, yscale = logx ? log10 : identity)
    xs = logx ? [r[:γ₂] for r in R] : dTs
    for (i, r) in enumerate(R); scatter!(a, [xs[i]], [y[i]], color = cols[i], marker = mk[i], markersize = 16, label = r[:label]); end
    if logx && pred !== nothing
        xa = [r[:γ₂] for r in A]; ya = y[[findfirst(==(r), R) for r in A]]
        xx = range(minimum(xa)*0.98, maximum(xa)*1.02, length = 20)
        lines!(a, xx, ya[1] .* (xx ./ xa[1]).^pred, color = :gray50, linestyle = :dash, label = @sprintf("inertial  γ₂^%+.2f", pred))
        fitexp !== nothing && lines!(a, xx, ya[1] .* (xx ./ xa[1]).^fitexp, color = :firebrick, label = @sprintf("A-series fit  γ₂^%+.2f", fitexp))
    end
    a
end
a = panel((1, 1), "mid-basin z20", "m", [-r[:z20m] for r in R]); axislegend(a, position = :rt, labelsize = 9)
panel((1, 2), "EUC core depth, 140W", "m", [-r[Symbol("140W")].zcore for r in R])
panel((1, 3), "EUC core speed, 140W (core band)", "m/s", [r[Symbol("140W")].core for r in R])
a = panel((2, 1), "mid-basin z20 vs γ₂", "m", [-r[:z20m] for r in R]; logx = true, pred = -0.5, fitexp = expH); axislegend(a, position = :rt, labelsize = 9)
panel((2, 2), "EUC core depth vs γ₂", "m", [-r[Symbol("140W")].zcore for r in R]; logx = true, pred = -0.5, fitexp = expHc)
panel((2, 3), "EUC core speed vs γ₂", "m/s", [r[Symbol("140W")].core for r in R]; logx = true, pred = 0.25, fitexp = expU)
save(joinpath(OUTDIR, "cmp_scaling.png"), fig); println("wrote cmp_scaling.png")
println("=== done → $OUTDIR")
