# =============================================================================
# plot_band_sections.jl — longitude–depth sections (latitude-band mean) of T, N² and u for
# the control and the warming branches, plus branch − control differences.
# Input: band_section_<N>deg.jld2 in each run folder (made by analysis/band_section.jl).
#
#     julia analysis/plot_band_sections.jl            (edit SET / DL / HALFDEG below)
# Writes <DL>/branch_comparison_<SET>/sections_<N>deg.png and sections_<N>deg_diff.png
# =============================================================================
using JLD2, Statistics, Printf, CairoMakie

DL  = expanduser("~/Downloads")
SET = length(ARGS) ≥ 1 ? ARGS[1] : "5y";  HALFDEG = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 5   # e.g. julia plot_band_sections.jl 1y 5
SETS = Dict(
    "1y" => [("control", "toy_model_control_2026-10-21_180403"), ("A1 +1 K", "toy_model_warmA1_2026-10-21_182550"),
             ("A2 +2 K", "toy_model_warmA2_2026-10-21_182936"), ("A4 +4 K", "toy_model_warmA4_2026-10-21_184307"),
             ("B2 +2 K ventilated", "toy_model_warmB2_2026-10-21_185704"), ("C2 +2/+4 K pattern", "toy_model_warmC2_2026-10-21_185857")],
    "5y" => [("control", "toy_model_control_2026-10-21_233303"), ("A1 +1 K", "toy_model_warmA1_2026-10-22_001758"),
             ("A2 +2 K", "toy_model_warmA2_2026-10-22_005915"), ("A4 +4 K", "toy_model_warmA4_2026-10-22_010512"),
             ("B2 +2 K ventilated", "toy_model_warmB2_2026-10-22_011949"), ("C2 +2/+4 K pattern", "toy_model_warmC2_2026-10-22_013131")])
OUTDIR = joinpath(DL, "branch_comparison_$SET"); mkpath(OUTDIR)
ZMAX = 400.0                                                     # plot depth (m)
lon(x) = (l = 120 + 160 * x / 1.5e7; l > 180 ? l - 360 : l)      # basin = 120°E → 80°W
lonlabel(l) = l < 0 ? @sprintf("%.0f°W", -l) : @sprintf("%.0f°E", l)

R = []
for (lab, folder) in SETS[SET]
    p = joinpath(DL, folder, @sprintf("band_section_%ddeg.jld2", HALFDEG))
    isfile(p) || error("missing $p — run analysis/band_section.jl on that folder first (where its toy_model.jld2 is)")
    d = load(p); push!(R, (; label = lab, x = d["x"], zC = d["zC"], zI = d["zI"], T = d["T"], N2 = d["N2"], u = d["u"], t0 = d["tstart"], t1 = d["tend"]))
end
C = R[1]; xl = lon.(C.x); xt = [lon(x) for x in range(0, 1.5e7, length = 9)]
ticks = (collect(range(0, 1.5e7, length = 9)) ./ 1e3, lonlabel.(xt))
xkm = C.x ./ 1e3
kT = C.zC .>= -ZMAX; kN = C.zI .>= -ZMAX

# ---------------- figure 1: fields ----------------
fig = Figure(size = (2000, 300 * length(R) + 120))
Label(fig[0, 1:6], @sprintf("|lat| ≤ %d° band mean, days %.0f–%.0f of each branch  —  T, N², u (longitude vs depth)", HALFDEG, C.t0, C.t1), fontsize = 20, tellwidth = false)
for (i, r) in enumerate(R)
    a = Axis(fig[i, 1], title = "$(r.label):  T (°C)", ylabel = "z (m)", xticks = ticks, xticklabelrotation = 0.5)
    h1 = heatmap!(a, xkm, r.zC[kT], r.T[:, kT], colormap = :thermal, colorrange = (10, 30)); contour!(a, xkm, r.zC[kT], r.T[:, kT], levels = [20.0], color = :white, linewidth = 1.5)
    a = Axis(fig[i, 3], title = "$(r.label):  N² (s⁻²)", xticks = ticks, xticklabelrotation = 0.5)
    h2 = heatmap!(a, xkm, r.zI[kN], r.N2[:, kN], colormap = :viridis, colorrange = (0, 4e-4))
    a = Axis(fig[i, 5], title = "$(r.label):  u (m/s)", xticks = ticks, xticklabelrotation = 0.5)
    # a ±5° band mean dilutes the EUC (half-width ~1.5°) about 3×, so the range is tighter than for core-band plots
    h3 = heatmap!(a, xkm, r.zC[kT], r.u[:, kT], colormap = :balance, colorrange = (-0.5, 0.5)); contour!(a, xkm, r.zC[kT], r.u[:, kT], levels = [0.0], color = :black, linewidth = 0.8)
    i == 1 && (Colorbar(fig[1, 2], h1, label = "°C"); Colorbar(fig[1, 4], h2, label = "s⁻²"); Colorbar(fig[1, 6], h3, label = "m/s"))
end
save(joinpath(OUTDIR, @sprintf("sections_%ddeg.png", HALFDEG)), fig); println("wrote sections_$(HALFDEG)deg.png")

# ---------------- figure 2: differences from control ----------------
fig = Figure(size = (2000, 300 * (length(R) - 1) + 120))
Label(fig[0, 1:6], @sprintf("|lat| ≤ %d° band mean: branch − control.  Black = control 20 °C isotherm and u = 0; colour = branch", HALFDEG), fontsize = 20, tellwidth = false)
for (i, r) in enumerate(R[2:end])
    a = Axis(fig[i, 1], title = "$(r.label):  ΔT (K)", ylabel = "z (m)", xticks = ticks, xticklabelrotation = 0.5)
    h1 = heatmap!(a, xkm, r.zC[kT], r.T[:, kT] .- C.T[:, kT], colormap = :balance, colorrange = (-4, 4))
    contour!(a, xkm, C.zC[kT], C.T[:, kT], levels = [20.0], color = :black, linewidth = 1.2); contour!(a, xkm, r.zC[kT], r.T[:, kT], levels = [20.0], color = :orangered, linewidth = 1.2)
    a = Axis(fig[i, 3], title = "$(r.label):  ΔN² (s⁻²)", xticks = ticks, xticklabelrotation = 0.5)
    h2 = heatmap!(a, xkm, r.zI[kN], r.N2[:, kN] .- C.N2[:, kN], colormap = :balance, colorrange = (-1.5e-4, 1.5e-4))
    a = Axis(fig[i, 5], title = "$(r.label):  Δu (m/s)", xticks = ticks, xticklabelrotation = 0.5)
    h3 = heatmap!(a, xkm, r.zC[kT], r.u[:, kT] .- C.u[:, kT], colormap = :balance, colorrange = (-0.2, 0.2))
    contour!(a, xkm, C.zC[kT], C.u[:, kT], levels = [0.0], color = :black, linewidth = 1.2); contour!(a, xkm, r.zC[kT], r.u[:, kT], levels = [0.0], color = :orangered, linewidth = 1.2)
    i == 1 && (Colorbar(fig[1, 2], h1, label = "K"); Colorbar(fig[1, 4], h2, label = "s⁻²"); Colorbar(fig[1, 6], h3, label = "m/s"))
end
save(joinpath(OUTDIR, @sprintf("sections_%ddeg_diff.png", HALFDEG)), fig); println("wrote sections_$(HALFDEG)deg_diff.png")
