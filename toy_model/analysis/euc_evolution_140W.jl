# Δu(z) and ΔT(z) at 140W, 5°S–5°N band mean, branch − control, by YEAR of the 5-year branches
# (the 2-day slice files hold the whole 5 years) — the toy-model analogue of OM2 Fig. 6a.
using JLD2, Statistics, Printf, CairoMakie
DL = expanduser("~/Downloads")
runs = [("control", "toy_model_control_2026-10-21_233303", 0.0), ("A1 +1 K", "toy_model_warmA1_2026-10-22_001758", 1.0), ("A2 +2 K", "toy_model_warmA2_2026-10-22_005915", 2.0),
        ("A4 +4 K", "toy_model_warmA4_2026-10-22_010512", 4.0), ("B2 +2 K ventilated", "toy_model_warmB2_2026-10-22_011949", 2.0), ("C2 +2/+4 K pattern", "toy_model_warmC2_2026-10-22_013131", 3.0)]
Nx, Ny, Nz = 600, 150, 60; Ly, Lz = 3.33e6, 2000.0; refinement, stretching = 3.5, 18; H = 4
hk(k) = (k-1)/Nz; ζ₀(k) = 1 + (hk(k)-1)/refinement; Σ(k) = (1-exp(-stretching*hk(k)))/(1-exp(-stretching))
zF = [Lz*(ζ₀(k)*Σ(k)-1) for k in 1:Nz+1]; zC = 0.5 .* (zF[1:end-1] .+ zF[2:end])
yC = [-Ly/2 + (j-0.5)*Ly/Ny for j in 1:Ny]; jb = findall(abs.(yC) .<= 5*111e3)
strip(A) = A[ntuple(d -> size(A, d) > 1 ? (H+1:size(A, d)-H) : (1:1), ndims(A))...]
yearly = Dict()
for (lab, folder, dT) in runs
    f = jldopen(joinpath(DL, folder, "toy_model_140W.jld2"), "r")
    its = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/u"])))); t = [f["timeseries/t/$k"] for k in its] ./ 86400 .- 1800
    U = zeros(5, Nz); T = zeros(5, Nz); n = zeros(Int, 5)
    for (i, k) in enumerate(its)
        y = clamp(ceil(Int, t[i]/365), 1, 5); t[i] <= 0 && continue
        U[y, :] .+= dropdims(mean(strip(f["timeseries/u/$k"])[1, jb, :]; dims = 1); dims = 1)
        T[y, :] .+= dropdims(mean(strip(f["timeseries/T/$k"])[1, jb, :]; dims = 1); dims = 1); n[y] += 1
    end
    close(f); yearly[lab] = (; u = U ./ n, T = T ./ n, dT)
    println("$lab: samples per year ", n)
end
C = yearly["control"]
println("\n140W, 5°S–5°N band: branch − control, yearly means.  (OM2 4° warming, Fig 6a: Δu ≈ +0.04 m/s at 90–100 m, −0.01 at 140–180 m, ≈ +0.01 at the surface after 36 yr of trend ≈ +1.9 K)")
println("branch           year   Δu(10m)  Δu(50m)  Δu(90m) Δu(140m) Δu(180m) | ΔT(10m) ΔT(50m) ΔT(100m) ΔT(150m)")
d(z) = argmin(abs.(zC .+ z))
for (lab, _, _) in runs[2:end], y in 1:5
    r = yearly[lab]; du = r.u[y, :] .- C.u[y, :]; dT = r.T[y, :] .- C.T[y, :]
    @printf "%-16s %4d   %+7.3f  %+7.3f  %+7.3f  %+7.3f  %+7.3f | %+6.2f  %+6.2f  %+6.2f  %+6.2f\n" lab y du[d(10)] du[d(50)] du[d(90)] du[d(140)] du[d(180)] dT[d(10)] dT[d(50)] dT[d(100)] dT[d(150)]
end
@printf "\ncontrol 140W 5°S–5°N mean u: surface %+.2f, core %+.2f at %.0f m (OM2 control: core ≈ +0.18 at ~120 m, surface ≈ −0.25)\n" C.u[5, end] maximum(C.u[5, :]) zC[argmax(C.u[5, :])]
fig = Figure(size = (1700, 650))
Label(fig[0, 1:5], "140W, 5°S–5°N band mean: Δu(z) and ΔT(z), branch − control, by year of the branch  (cf. OM2 Fig. 6a: trend peak at 90–100 m, negative below ~130 m)", fontsize = 16, tellwidth = false)
ycols = [:gold, :orange, :orangered, :crimson, :darkred]
for (i, (lab, _, _)) in enumerate(runs[2:end])
    r = yearly[lab]
    a = Axis(fig[1, i], title = lab, xlabel = "Δu (m/s)", ylabel = i == 1 ? "z (m)" : "")
    for y in 1:5; lines!(a, r.u[y, :] .- C.u[y, :], zC, color = ycols[y], label = "year $y"); end
    vlines!(a, [0.0], color = :black, linestyle = :dash); ylims!(a, -350, 0); xlims!(a, -0.2, 0.15)
    hlines!(a, [zC[argmax(C.u[5, :])]], color = :gray, linestyle = :dot)
    i == 1 && axislegend(a, position = :lb, labelsize = 9)
end
save(joinpath(DL, "branch_comparison_5y", "euc_evolution_140W_5deg.png"), fig); println("wrote euc_evolution_140W_5deg.png")
