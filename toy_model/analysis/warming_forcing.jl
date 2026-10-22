# =============================================================================
# warming_forcing.jl — the surface warming applied to the model: the SST restoring target
# along the basin for the control and each warming branch.
#
#     julia analysis/warming_forcing.jl            (edit RUNS / DL below)
# Writes <DL>/branch_comparison_5y/warming_forcing.png
# Reads toy_model_surface.jld2 (2-day top-cell slices); last 365 days.
# =============================================================================
using JLD2, Statistics, Printf, CairoMakie

DL   = expanduser("~/Downloads")
RUNS = [("control", "toy_model_control_2026-10-21_233303", 0.0, :black),
        ("+1 K",    "toy_model_warmA1_2026-10-22_001758", 1.0, :dodgerblue),
        ("+2 K",    "toy_model_warmA2_2026-10-22_005915", 2.0, :darkorange),
        ("+4 K",    "toy_model_warmA4_2026-10-22_010512", 4.0, :crimson)]
OUT  = joinpath(DL, "branch_comparison_5y"); mkpath(OUT)
WINDOW = 365.0

# restoring (toy_model.jl): F = λ_T (T − T_target(x)),  T_target = T_cold + (T_warm − T_cold)·½[1 − tanh((x/Lx − x_f)/w_f)]
const Nx, Ny = 600, 150;  const Lx, Ly = 1.5e7, 3.33e6;  const H = 4
λ_T = 4e-5; T_warm0 = 29.5; T_cold0 = 24.0; x_f = 0.6; w_f = 0.18; ρ₀ = 1026.0; cₚ = 3991.0
xC = [(i-0.5)*Lx/Nx for i in 1:Nx];  yC = [-Ly/2 + (j-0.5)*Ly/Ny for j in 1:Ny]
lon(x) = (l = 120 + 160x/Lx; l > 180 ? l - 360 : l);  lonlabel(l) = l < 0 ? @sprintf("%.0f°W", -l) : @sprintf("%.0f°E", l)
xticks = (collect(range(0, Lx, length = 9)) ./ 1e3, lonlabel.(lon.(range(0, Lx, length = 9))))
Ttarget(x, ΔT) = (T_cold0 + ΔT) + (T_warm0 - T_cold0) * 0.5 * (1 - tanh((x/Lx - x_f)/w_f))
strip(A) = A[ntuple(d -> size(A, d) > 1 ? (H+1:size(A, d)-H) : (1:1), ndims(A))...]
jeq = findall(abs.(yC) .<= 111e3)      # ±1°

R = []
for (lab, folder, ΔT, col) in RUNS
    f = jldopen(joinpath(DL, folder, "toy_model_surface.jld2"), "r")
    its = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/T"])))); t = [f["timeseries/t/$k"] for k in its] ./ 86400
    sel = its[t .>= t[end] - WINDOW]
    SST = mean([strip(f["timeseries/T/$k"])[:, :, 1] for k in sel]); close(f)
    tar = Ttarget.(xC, ΔT)
    push!(R, (; lab, ΔT, col, SST, sst_eq = dropdims(mean(SST[:, jeq]; dims = 2); dims = 2), tar,
              Q_eq = ρ₀*cₚ*λ_T .* (tar .- dropdims(mean(SST[:, jeq]; dims = 2); dims = 2))))   # W/m², + into the ocean
end
C = R[1]
for r in R; @printf "%-8s target W/E %.1f/%.1f | SST eq W/E %.1f/%.1f | basin-mean ΔSST vs control %+.2f K\n" r.lab r.tar[1] r.tar[end] r.sst_eq[1] r.sst_eq[end] mean(r.SST .- C.SST); end

fig = Figure(size = (900, 500), fontsize = 16)
a = Axis(fig[1, 1], title = "Surface warming: SST restoring target", ylabel = "°C", xticks = xticks)
for r in R; lines!(a, xC ./ 1e3, r.tar, color = r.col, linewidth = 3, label = r.lab); end
axislegend(a, position = :lb)
save(joinpath(OUT, "warming_forcing.png"), fig); println("wrote ", joinpath(OUT, "warming_forcing.png"))
