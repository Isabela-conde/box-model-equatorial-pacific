# =============================================================================
# warming_profiles_140W.jl — surface-warming branches vs control at 140°W, 5°S–5°N mean:
# u(z), Δu(z) and ΔT(z).  Dashed line = thermocline (depth of max ∂b/∂z of the band-mean
# profile), dot = core (max u), both in the run's colour.
#
#     julia analysis/warming_profiles_140W.jl          (edit RUNS / DL / BAND below)
# Writes <DL>/branch_comparison_5y/warming_profiles_140W.png
# Reads toy_model_140W.jld2 (2-day y–z slices) from each run folder; last 365 days.
# =============================================================================
using JLD2, Statistics, Printf, CairoMakie, SeawaterPolynomials, SeawaterPolynomials.TEOS10

DL   = expanduser("~/Downloads")
RUNS = [("control", "toy_model_control_2026-10-21_233303", :black),
        ("+1 K",    "toy_model_warmA1_2026-10-22_001758", :dodgerblue),
        ("+2 K",    "toy_model_warmA2_2026-10-22_005915", :darkorange),
        ("+4 K",    "toy_model_warmA4_2026-10-22_010512", :crimson)]
OUT  = joinpath(DL, "branch_comparison_5y"); mkpath(OUT)
BAND = 5.0          # half-width of the latitude band (degrees)
WINDOW = 365.0

const Nx, Ny, Nz = 600, 150, 60;  const Ly, Lz = 3.33e6, 2000.0
const refinement, stretching = 3.5, 18;  const H = 4;  const g, ρ₀ = 9.81, 1026.0
hk(k) = (k-1)/Nz;  ζ₀(k) = 1 + (hk(k)-1)/refinement;  Σ(k) = (1 - exp(-stretching*hk(k)))/(1 - exp(-stretching))
zF = [Lz*(ζ₀(k)*Σ(k) - 1) for k in 1:Nz+1];  zC = 0.5 .* (zF[1:end-1] .+ zF[2:end]);  zI = zF[2:end-1]
yC = [-Ly/2 + (j-0.5)*Ly/Ny for j in 1:Ny];  jb = findall(abs.(yC) .<= BAND*111e3)
strip(A) = A[ntuple(d -> size(A, d) > 1 ? (H+1:size(A, d)-H) : (1:1), ndims(A))...]
eos = TEOS10EquationOfState(reference_density = ρ₀)
thermocline(T, S) = zI[argmax(diff([-g*SeawaterPolynomials.ρ′(T[k], S[k], zC[k], eos)/ρ₀ for k in 1:Nz]) ./ diff(zC))]
kz = findall(-350 .<= zC .<= -5)
core(u) = (k = kz[argmax(u[kz])]; (u[k], zC[k]))

R = []
for (lab, folder, col) in RUNS
    f = jldopen(joinpath(DL, folder, "toy_model_140W.jld2"), "r")
    its = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/u"])))); t = [f["timeseries/t/$k"] for k in its] ./ 86400
    sel = its[t .>= t[end] - WINDOW]
    p(v) = mean([dropdims(mean(strip(f["timeseries/$v/$k"])[1, jb, :]; dims = 1); dims = 1) for k in sel])
    u = p("u"); T = p("T"); S = p("S"); close(f)
    push!(R, (; lab, col, u, T, zt = thermocline(T, S)))
end
C = R[1]
println("thermocline depth: ", [(r.lab, round(r.zt)) for r in R])

bandlab = @sprintf("%.0f°S–%.0f°N mean", BAND, BAND)
fig = Figure(size = (1800, 620), fontsize = 16)
Label(fig[0, 1:3], "Profiles at 140°W, $bandlab", fontsize = 20, tellwidth = false)
a = Axis(fig[1, 1], title = "u(z)", xlabel = "u (m/s)", ylabel = "depth (m)")
for r in R
    lines!(a, r.u, zC, color = r.col, linewidth = 3, label = r.lab)
    hlines!(a, [r.zt], color = r.col, linestyle = :dash, linewidth = 1.5)
    uc, zc = core(r.u); scatter!(a, [uc], [zc], color = r.col, markersize = 16, strokecolor = :white, strokewidth = 1)
end
vlines!(a, [0.0], color = :gray50); ylims!(a, -350, 0); axislegend(a, position = :lt)
a = Axis(fig[1, 2], title = "Δu(z)", xlabel = "Δu (m/s)")
for r in R; hlines!(a, [r.zt], color = r.col, linestyle = :dash, linewidth = 1.5); end
for r in R[2:end]; lines!(a, r.u .- C.u, zC, color = r.col, linewidth = 3, label = r.lab); end
vlines!(a, [0.0], color = :gray50); ylims!(a, -350, 0); axislegend(a, position = :rb)
a = Axis(fig[1, 3], title = "ΔT(z)", xlabel = "ΔT (K)")
for r in R; hlines!(a, [r.zt], color = r.col, linestyle = :dash, linewidth = 1.5); end
for r in R[2:end]; lines!(a, r.T .- C.T, zC, color = r.col, linewidth = 3, label = r.lab); end
vlines!(a, [0.0], color = :gray50); ylims!(a, -350, 0); axislegend(a, position = :rb)
save(joinpath(OUT, "warming_profiles_140W.png"), fig); println("wrote ", joinpath(OUT, "warming_profiles_140W.png"))
