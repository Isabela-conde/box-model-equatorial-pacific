# =============================================================================
# warming_yz_140W.jl — surface-warming branches vs control at 140°W, latitude–depth:
# control u(y,z) and Δu(y,z) for each warming amplitude, with the thermocline
# (depth of max ∂b/∂z) overlaid: black = control, red = warm run.
#
#     julia analysis/warming_yz_140W.jl          (edit RUNS / DL below)
# Writes <DL>/branch_comparison_5y/warming_yz_140W.png
# Reads toy_model_140W.jld2 (2-day y–z slices) from each run folder; last 365 days.
# =============================================================================
using JLD2, Statistics, Printf, CairoMakie, SeawaterPolynomials, SeawaterPolynomials.TEOS10

DL   = expanduser("~/Downloads")
RUNS = [("control", "toy_model_control_2026-10-21_233303"),
        ("+1 K",    "toy_model_warmA1_2026-10-22_001758"),
        ("+2 K",    "toy_model_warmA2_2026-10-22_005915"),
        ("+4 K",    "toy_model_warmA4_2026-10-22_010512")]
OUT  = joinpath(DL, "branch_comparison_5y"); mkpath(OUT)
WINDOW = 365.0

const Nx, Ny, Nz = 600, 150, 60;  const Ly, Lz = 3.33e6, 2000.0
const refinement, stretching = 3.5, 18;  const H = 4;  const g, ρ₀ = 9.81, 1026.0
hk(k) = (k-1)/Nz;  ζ₀(k) = 1 + (hk(k)-1)/refinement;  Σ(k) = (1 - exp(-stretching*hk(k)))/(1 - exp(-stretching))
zF = [Lz*(ζ₀(k)*Σ(k) - 1) for k in 1:Nz+1];  zC = 0.5 .* (zF[1:end-1] .+ zF[2:end]);  zI = zF[2:end-1]
yC = [-Ly/2 + (j-0.5)*Ly/Ny for j in 1:Ny];  lat = yC ./ 111e3
strip(A) = A[ntuple(d -> size(A, d) > 1 ? (H+1:size(A, d)-H) : (1:1), ndims(A))...]
eos = TEOS10EquationOfState(reference_density = ρ₀)

# thermocline per latitude = depth of max ∂b/∂z, lightly smoothed in latitude
function thermocline(T, S)
    zt = [zI[argmax(diff([-g*SeawaterPolynomials.ρ′(T[j, k], S[j, k], zC[k], eos)/ρ₀ for k in 1:Nz]) ./ diff(zC))] for j in 1:Ny]
    [mean(zt[max(1, j-1):min(Ny, j+1)]) for j in 1:Ny]
end

R = []
for (lab, folder) in RUNS
    f = jldopen(joinpath(DL, folder, "toy_model_140W.jld2"), "r")
    its = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/u"])))); t = [f["timeseries/t/$k"] for k in its] ./ 86400
    sel = its[t .>= t[end] - WINDOW]
    U = mean([strip(f["timeseries/u/$k"])[1, :, :] for k in sel]); T = mean([strip(f["timeseries/T/$k"])[1, :, :] for k in sel]); S = mean([strip(f["timeseries/S/$k"])[1, :, :] for k in sel])
    close(f); push!(R, (; lab, U, zt = thermocline(T, S)))
end
C = R[1]

fig = Figure(size = (2000, 640), fontsize = 17)
Label(fig[0, 1:4], "Surface warming at 140°W: zonal velocity across latitude", fontsize = 18, tellwidth = false)
zlim = (-350, 0); ylim = (-6, 6)
a = Axis(fig[1, 1], title = "control  u (m/s)", xlabel = "latitude (°)", ylabel = "depth (m)")
hm = heatmap!(a, lat, zC, C.U, colormap = :balance, colorrange = (-1.0, 1.0)); contour!(a, lat, zC, C.U, levels = 0.2:0.2:1.0, color = :black, linewidth = 0.7)
lines!(a, lat, C.zt, color = :black, linestyle = :dash, linewidth = 2.5)
xlims!(a, ylim...); ylims!(a, zlim...)
Colorbar(fig[2, 1], hm, label = "u (m/s)", vertical = false, flipaxis = false)
hm2 = nothing
for (i, r) in enumerate(R[2:end])
    a = Axis(fig[1, i + 1], title = "$(r.lab):  Δu (m/s)", xlabel = "latitude (°)")
    global hm2 = heatmap!(a, lat, zC, r.U .- C.U, colormap = :balance, colorrange = (-0.15, 0.15)); contour!(a, lat, zC, C.U, levels = 0.2:0.2:1.0, color = :black, linewidth = 0.7)
    lines!(a, lat, C.zt, color = :black, linestyle = :dash, linewidth = 2.5); lines!(a, lat, r.zt, color = :red, linestyle = :dash, linewidth = 2.5)
    xlims!(a, ylim...); ylims!(a, zlim...); hideydecorations!(a, grid = false)
end
Colorbar(fig[2, 2:4], hm2, label = "Δu (m/s)", vertical = false, flipaxis = false, width = Relative(0.6))
save(joinpath(OUT, "warming_yz_140W.png"), fig); println("wrote ", joinpath(OUT, "warming_yz_140W.png"))
