# Equatorial zonal pressure-gradient force from the two saved y–z sections (140W and 110W, full T,S):
# p(y,z) = g∫_z^0 ρ' dz' + ρ₀ g η_s with η_s the steric height referenced to 2000 m (level of no
# motion).  PGF(y,z) = −(p_110W − p_140W)/(ρ₀ Δx) is the force per unit mass averaged over 140–110W.
# Compared with Δu(y,z) at 140W.  Last year of the 5-yr branches.
using JLD2, Statistics, Printf, CairoMakie, SeawaterPolynomials, SeawaterPolynomials.TEOS10
DL = expanduser("~/Downloads")
runs = [("control", "toy_model_control_2026-10-21_233303", :black), ("+1 K", "toy_model_warmA1_2026-10-22_001758", :dodgerblue),
        ("+2 K", "toy_model_warmA2_2026-10-22_005915", :darkorange), ("+4 K", "toy_model_warmA4_2026-10-22_010512", :crimson)]
Nx, Ny, Nz = 600, 150, 60; Lx, Ly, Lz = 1.5e7, 3.33e6, 2000.0; refinement, stretching = 3.5, 18; H = 4
hk(k) = (k-1)/Nz; ζ₀(k) = 1 + (hk(k)-1)/refinement; Σ(k) = (1-exp(-stretching*hk(k)))/(1-exp(-stretching))
zF = [Lz*(ζ₀(k)*Σ(k)-1) for k in 1:Nz+1]; zC = 0.5 .* (zF[1:end-1] .+ zF[2:end]); dz = diff(zF)
yC = [-Ly/2 + (j-0.5)*Ly/Ny for j in 1:Ny]; lat = yC ./ 111e3
strip(A) = A[ntuple(d -> size(A, d) > 1 ? (H+1:size(A, d)-H) : (1:1), ndims(A))...]
eos = TEOS10EquationOfState(reference_density = 1026); g = 9.81; ρ₀ = 1026.0
Δx = (0.81 - 0.625) * Lx
function pressure(T, S)                                   # (Ny, Nz) → p anomaly (Pa) incl. steric surface
    ρ = [SeawaterPolynomials.ρ′(Float64(T[j, k]), Float64(S[j, k]), zC[k], eos) for j in 1:Ny, k in 1:Nz]
    η = [-(1/ρ₀)*sum(ρ[j, :] .* dz) for j in 1:Ny]
    [g*sum(ρ[j, k:Nz] .* dz[k:Nz]) + ρ₀*g*η[j] for j in 1:Ny, k in 1:Nz]
end
R = []
for (lab, folder, col) in runs
    P = Dict(); U = Dict()
    for sec in ("140W", "110W")
        f = jldopen(joinpath(DL, folder, "toy_model_$sec.jld2"), "r")
        its = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/u"])))); t = [f["timeseries/t/$k"] for k in its] ./ 86400; sel = its[t .>= t[end] - 365]
        T = mean([strip(f["timeseries/T/$k"])[1, :, :] for k in sel]); S = mean([strip(f["timeseries/S/$k"])[1, :, :] for k in sel]); U[sec] = mean([strip(f["timeseries/u/$k"])[1, :, :] for k in sel]); close(f)
        P[sec] = pressure(T, S)
    end
    pgf = -(P["110W"] .- P["140W"]) ./ (ρ₀ * Δx)
    push!(R, (; lab, col, pgf, u = U["140W"]))
end
C = R[1]; jeq = findall(abs.(yC) .<= 100e3); d(z) = argmin(abs.(zC .+ z))
println("Equator (|y| ≤ 100 km), PGF averaged 140–110W, 10⁻⁷ m s⁻² ; u at 140W")
println("run       PGF 25m   50m   75m  100m  150m  200m | u(140W) 25m   50m  100m  150m")
for r in R
    pe = dropdims(mean(r.pgf[jeq, :]; dims = 1); dims = 1); ue = dropdims(mean(r.u[jeq, :]; dims = 1); dims = 1)
    @printf "%-8s   %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f | %+5.2f %+5.2f %+5.2f %+5.2f\n" r.lab 1e7*pe[d(25)] 1e7*pe[d(50)] 1e7*pe[d(75)] 1e7*pe[d(100)] 1e7*pe[d(150)] 1e7*pe[d(200)] ue[d(25)] ue[d(50)] ue[d(100)] ue[d(150)]
end
println("Δ(warm − control):  ΔPGF (%)                              | Δu (m/s)")
for r in R[2:end]
    pe = dropdims(mean(r.pgf[jeq, :]; dims = 1); dims = 1); pc = dropdims(mean(C.pgf[jeq, :]; dims = 1); dims = 1)
    ue = dropdims(mean(r.u[jeq, :]; dims = 1); dims = 1); uc = dropdims(mean(C.u[jeq, :]; dims = 1); dims = 1)
    @printf "%-8s   %+5.0f %+5.0f %+5.0f %+5.0f %+5.0f %+5.0f | %+5.2f %+5.2f %+5.2f %+5.2f\n" r.lab [100*(pe[d(z)]/pc[d(z)]-1) for z in (25,50,75,100,150,200)]... ue[d(25)]-uc[d(25)] ue[d(50)]-uc[d(50)] ue[d(100)]-uc[d(100)] ue[d(150)]-uc[d(150)]
end
fig = Figure(size = (2000, 900), fontsize = 16)
Label(fig[0, 1:4], "Zonal pressure-gradient force (140–110°W mean) and Δu at 140°W, last year", fontsize = 18, tellwidth = false)
a = Axis(fig[1, 1], title = "control PGF (10⁻⁷ m s⁻²)", ylabel = "depth (m)"); hm = heatmap!(a, lat, zC, 1e7 .* C.pgf, colormap = :balance, colorrange = (-6, 6)); contour!(a, lat, zC, C.u, levels = 0.2:0.2:1.0, color = :black, linewidth = 0.7); xlims!(a, -6, 6); ylims!(a, -350, 0); Colorbar(fig[1, 0], hm, flipaxis = false)
for (i, r) in enumerate(R[2:end])
    a = Axis(fig[1, i+1], title = "$(r.lab): ΔPGF (10⁻⁷ m s⁻²)"); hm2 = heatmap!(a, lat, zC, 1e7 .* (r.pgf .- C.pgf), colormap = :balance, colorrange = (-1.5, 1.5)); contour!(a, lat, zC, C.u, levels = 0.2:0.2:1.0, color = :black, linewidth = 0.7); xlims!(a, -6, 6); ylims!(a, -350, 0); hideydecorations!(a, grid = false)
    i == 3 && Colorbar(fig[1, 5], hm2)
    a = Axis(fig[2, i+1], title = "$(r.lab): Δu at 140°W (m/s)", xlabel = "latitude (°)"); hm3 = heatmap!(a, lat, zC, r.u .- C.u, colormap = :balance, colorrange = (-0.15, 0.15)); contour!(a, lat, zC, C.u, levels = 0.2:0.2:1.0, color = :black, linewidth = 0.7); xlims!(a, -6, 6); ylims!(a, -350, 0); hideydecorations!(a, grid = false)
    i == 3 && Colorbar(fig[2, 5], hm3)
end
a = Axis(fig[2, 1], title = "equator (|y| ≤ 100 km): ΔPGF (solid) and Δu (dashed)", xlabel = "ΔPGF (10⁻⁷ m s⁻²)", ylabel = "depth (m)")
a2 = Axis(fig[2, 1], xaxisposition = :top, xlabel = "Δu (m/s)"); hidespines!(a2); hideydecorations!(a2)
for r in R[2:end]
    lines!(a, 1e7 .* dropdims(mean(r.pgf[jeq, :] .- C.pgf[jeq, :]; dims = 1); dims = 1), zC, color = r.col, linewidth = 3)
    lines!(a2, dropdims(mean(r.u[jeq, :] .- C.u[jeq, :]; dims = 1); dims = 1), zC, color = r.col, linewidth = 2, linestyle = :dash)
end
vlines!(a, [0.0], color = :gray50); ylims!(a, -350, 0); ylims!(a2, -350, 0); xlims!(a, -1.5, 1.5); xlims!(a2, -0.3, 0.3)
save(joinpath(DL, "branch_comparison_5y", "pgf_sections.png"), fig); println("wrote pgf_sections.png")
