# Paper-1 §1/§3 metrics from the 5-yr branches (140W & 110W slices, last 365 d):
# centroid depth, peak speed, half-thickness of the eastward jet; core-relative N², S², Ri and the
# Ri-scheme κ implied by them (κ = κ₀·½[1−tanh((Ri−0.1)/0.4)] + κ_bg) at fixed distances above/below the core.
using JLD2, Statistics, Printf, SeawaterPolynomials, SeawaterPolynomials.TEOS10
DL = expanduser("~/Downloads")
runs = [("control", "toy_model_control_2026-10-21_233303", 0.0), ("A1 +1K", "toy_model_warmA1_2026-10-22_001758", 1.0), ("A2 +2K", "toy_model_warmA2_2026-10-22_005915", 2.0),
        ("A4 +4K", "toy_model_warmA4_2026-10-22_010512", 4.0), ("B2 vent", "toy_model_warmB2_2026-10-22_011949", 2.0), ("C2 patt", "toy_model_warmC2_2026-10-22_013131", 3.0)]
Nx, Ny, Nz = 600, 150, 60; Ly, Lz = 3.33e6, 2000.0; refinement, stretching = 3.5, 18; H = 4
hk(k) = (k-1)/Nz; ζ₀(k) = 1 + (hk(k)-1)/refinement; Σ(k) = (1-exp(-stretching*hk(k)))/(1-exp(-stretching))
zF = [Lz*(ζ₀(k)*Σ(k)-1) for k in 1:Nz+1]; zC = 0.5 .* (zF[1:end-1] .+ zF[2:end]); zI = zF[2:end-1]; dz = diff(zF)
yC = [-Ly/2 + (j-0.5)*Ly/Ny for j in 1:Ny]
strip(A) = A[ntuple(d -> size(A, d) > 1 ? (H+1:size(A, d)-H) : (1:1), ndims(A))...]
eos = TEOS10EquationOfState(reference_density = 1026); g = 9.81
taper(Ri) = 0.5*(1 - tanh((Ri - 0.1)/0.4)); κofRi(Ri) = 1e-2*taper(Ri) + 1e-5
interp(prof, zz, z0) = begin k = argmin(abs.(zz .- z0)); prof[k] end
for (band, halfkm) in (("core |y|≤100 km", 100e3), ("5°S–5°N", 555e3))
    jb = findall(abs.(yC) .<= halfkm)
    println("\n=== band $band, last 365 d.  centroid = ∫z·u⁺dz/∫u⁺dz over 0–350 m; thickness = depth range with u > ½·peak")
    println(@sprintf("%-8s %-5s %6s %7s %7s %7s | %7s %7s %7s %7s | %7s %7s %7s", "run", "sec", "peak", "z_peak", "z_cent", "thick", "N²c-30", "S²c-30", "Ri c-30", "κ c-30", "N²c+30", "S²c+30", "Ri c+30"))
    for sec in ("140W", "110W"), (lab, folder, dT) in runs
        f = jldopen(joinpath(DL, folder, "toy_model_$sec.jld2"), "r")
        its = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/u"])))); t = [f["timeseries/t/$k"] for k in its] ./ 86400; sel = its[t .>= t[end] - 365]
        p(v) = mean([dropdims(mean(strip(f["timeseries/$v/$k"])[1, jb, :]; dims = 1); dims = 1) for k in sel])
        u = p("u"); v = p("v"); T = p("T"); S = p("S"); close(f)
        kz = findall(-350 .<= zC .<= -5); up = max.(u[kz], 0.0)
        zcent = sum(zC[kz] .* up .* dz[kz]) / sum(up .* dz[kz]); kpk = kz[argmax(u[kz])]; U = u[kpk]
        half = findall(k -> u[k] > 0.5U, kz); thick = sum(dz[kz][half])
        ρ = [SeawaterPolynomials.ρ′(T[k], S[k], zC[k], eos) for k in 1:Nz]
        N2 = -(g/1026) .* diff(ρ) ./ diff(zC); S2 = (diff(u) ./ diff(zC)).^2 .+ (diff(v) ./ diff(zC)).^2; Ri = N2 ./ max.(S2, 1e-12)
        zc = zC[kpk]
        @printf "%-8s %-5s %6.2f %7.0f %7.0f %7.0f | %7.1e %7.1e %7.2f %7.1e | %7.1e %7.1e %7.2f\n" lab sec U zc zcent thick interp(N2, zI, zc+30) interp(S2, zI, zc+30) interp(Ri, zI, zc+30) κofRi(interp(Ri, zI, zc+30)) interp(N2, zI, zc-30) interp(S2, zI, zc-30) interp(Ri, zI, zc-30)
    end
end
