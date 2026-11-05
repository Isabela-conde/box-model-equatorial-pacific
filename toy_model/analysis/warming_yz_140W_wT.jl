# =============================================================================
# warming_yz_140W_wT.jl — same layout as warming_yz_140W.jl, for vertical velocity and
# temperature: control w(y,z) / T(y,z) and Δw / ΔT for each warming amplitude.
# Solid labelled contours = control mean u on every panel; thermocline (depth of max ∂b/∂z)
# overlaid: black = control, red = warm run.
#
#     julia analysis/warming_yz_140W_wT.jl          (edit RUNS / DL below)
# Writes <DL>/branch_comparison_5y/warming_yz_140W_w.png and warming_yz_140W_T.png
# Reads toy_model_140W.jld2 (2-day y–z slices) from each run folder; last 365 days.
# w is converted to m/day.
# =============================================================================
using JLD2, Statistics, Printf, CairoMakie, SeawaterPolynomials, SeawaterPolynomials.TEOS10

DL   = expanduser("~/Downloads")
RUNS = [("control", "toy_model_control_2026-10-21_233303"),
        ("+1 K",    "toy_model_warmA1_2026-10-22_001758"),
        ("+2 K",    "toy_model_warmA2_2026-10-22_005915"),
        ("+4 K",    "toy_model_warmA4_2026-10-22_010512")]
OUT  = joinpath(DL, "branch_comparison_5y"); mkpath(OUT)
WINDOW = 365.0

# colour ranges
W_RANGE  = (-3.0, 3.0)     # control w (m/day)
DW_RANGE = (-1.0, 1.0)     # Δw (m/day)
T_RANGE  = (8.0, 30.0)     # control T (°C)
DT_RANGE = (-4.0, 4.0)     # ΔT (K)

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
    U = mean([strip(f["timeseries/u/$k"])[1, :, :] for k in sel])
    W = mean([strip(f["timeseries/w/$k"])[1, :, :] for k in sel]) .* 86400   # m/day, on z faces (Nz+1)
    T = mean([strip(f["timeseries/T/$k"])[1, :, :] for k in sel])
    S = mean([strip(f["timeseries/S/$k"])[1, :, :] for k in sel])
    close(f); push!(R, (; lab, U, W, T, zt = thermocline(T, S)))
    @printf("%-8s  w: min %.2f max %.2f m/day   T: min %.1f max %.1f °C\n", lab, minimum(W), maximum(W), minimum(T), maximum(T))
end
C = R[1]

# control mean-u contours: solid for eastward, dashed for westward, labelled inline
ucontours!(ax, U) = begin
    contour!(ax, lat, zC, U, levels = [0.2, 0.4, 0.6, 0.8, 1.0], color = :black, linewidth = 1.2, labels = true, labelsize = 13, labelcolor = :black)
    contour!(ax, lat, zC, U, levels = [-0.4, -0.2], color = :black, linewidth = 1.2, linestyle = :dash, labels = true, labelsize = 13, labelcolor = :black)
end

# one figure: control field + Δfield per warming run.  `sym` = :W or :T, `z` = its vertical coordinate.
function make_figure(sym, z, name, units, suptitle, cmap, crange, dcmap, dcrange, outfile)
    fig = Figure(size = (1900, 550), fontsize = 17)
    Label(fig[0, 1:4], suptitle, fontsize = 18, tellwidth = false)
    zlim = (-350, 0); ylim = (-10, 10)
    a = Axis(fig[1, 1], title = "control  $name ($units)", xlabel = "latitude (°)", ylabel = "depth (m)")
    hm = heatmap!(a, lat, z, getfield(C, sym), colormap = cmap, colorrange = crange)
    ucontours!(a, C.U)
    lines!(a, lat, C.zt, color = :black, linestyle = :dash, linewidth = 2.5)
    xlims!(a, ylim...); ylims!(a, zlim...)
    Colorbar(fig[2, 1], hm, label = "$name ($units)", vertical = false, flipaxis = false)
    hm2 = nothing
    for (i, r) in enumerate(R[2:end])
        a = Axis(fig[1, i + 1], title = "$(r.lab):  Δ$name ($units)", xlabel = "latitude (°)")
        hm2 = heatmap!(a, lat, z, getfield(r, sym) .- getfield(C, sym), colormap = dcmap, colorrange = dcrange)
        ucontours!(a, C.U)   # control mean u
        lines!(a, lat, C.zt, color = :black, linestyle = :dash, linewidth = 2.5); lines!(a, lat, r.zt, color = :red, linestyle = :dash, linewidth = 2.5)
        xlims!(a, ylim...); ylims!(a, zlim...); hideydecorations!(a, grid = false)
    end
    Colorbar(fig[2, 2:4], hm2, label = "Δ$name ($units)", vertical = false, flipaxis = false, width = Relative(0.6))
    save(joinpath(OUT, outfile), fig); println("wrote ", joinpath(OUT, outfile))
end

make_figure(:W, zF, "w", "m/day", "Surface warming at 140°W: vertical velocity across latitude",
            :balance, W_RANGE, :balance, DW_RANGE, "warming_yz_140W_w.png")
make_figure(:T, zC, "T", "°C",    "Surface warming at 140°W: temperature across latitude",
            :thermal, T_RANGE, :balance, DT_RANGE, "warming_yz_140W_T.png")
