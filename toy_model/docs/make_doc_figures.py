"""Generate boundary-condition / forcing figures for the toy-model description.
All curves are computed directly from the model's functional forms in toy_model.jl
(production configuration of 2026-10-21).  Figures are written to docs/figures/.

Run from anywhere:  python3 toy_model/docs/make_doc_figures.py
Figure 08 (result section) is drawn from the production spin-up if its 3-D output
file can be found (see RESULT_RUNS below); otherwise it is skipped."""
import numpy as np
import matplotlib.pyplot as plt
import os, h5py

HERE = os.path.dirname(os.path.abspath(__file__))
FIG  = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

# ---- domain / constants (from toy_model.jl) ----
Lx, Ly, Lz = 1.5e7, 3.33e6, 2000.0           # ±15° domain
Nx, Ny, Nz = 600, 150, 60                    # Δx ≈ 25 km, Δy ≈ 22 km
rho0, cp = 1026.0, 3991.0
deg = 111e3                                  # m per degree latitude

def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, name), dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)

# ================= 1. vertical grid =================
ref, st = 3.5, 18
def zif(k):
    hh = (k-1)/Nz; z0 = 1+(hh-1)/ref
    S = (1-np.exp(-st*hh))/(1-np.exp(-st)); return Lz*(z0*S-1)
zf = np.array([zif(k) for k in range(1, Nz+2)])
zc = 0.5*(zf[:-1]+zf[1:]); dz = np.diff(zf)
fig, ax = plt.subplots(1, 2, figsize=(9, 5))
ax[0].plot(dz, zc, "-o", ms=3); ax[0].set_xlabel("Δz (m)"); ax[0].set_ylabel("z (m)")
ax[0].set_title("Vertical cell thickness"); ax[0].grid(alpha=.3)
ax[1].plot(dz, zc, "-o", ms=3); ax[1].set_ylim(-400, 0); ax[1].set_xlabel("Δz (m)")
ax[1].set_title("Upper 400 m (≈9.5 m uniform)"); ax[1].grid(alpha=.3)
save(fig, "01_vertical_grid.png")

# ================= 2. wind stress =================
# Matches surface_u10_profile in toy_model.jl: equator dip u10_eq, off-equatorial
# bump (lat/jet_lat)²·exp(1-(lat/jet_lat)²) toward u10_trades at ±jet_lat, minus a
# (lat/lat_edge)⁴ wall taper that weakens the easterlies to a NON-ZERO value at the walls.
u10_eq, u10_trades, jet_lat, wall_drop, cD, rhoa = 7.0, 8.0, 10.0, 4.0, 1.2e-3, 1.2
lat_edge = (Ly/2)/deg
y = np.linspace(-Ly/2, Ly/2, 600); lat = y/deg
trade = (lat/jet_lat)**2 * np.exp(1 - (lat/jet_lat)**2)
edge  = (lat/lat_edge)**4
speed = u10_eq + (u10_trades - u10_eq)*trade - wall_drop*edge   # trade-wind SPEED
u10   = -speed                               # easterly trades → NEGATIVE zonal wind
tau_x = -rhoa*cD*speed**2                    # physical eastward stress → NEGATIVE (westward)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(u10, lat); ax[0].set_xlabel(r"zonal wind $u_{10}$ (m s$^{-1}$)   [$-$ = easterly]")
ax[0].set_ylabel("latitude (°)"); ax[0].set_title("Surface zonal wind (easterly trades)")
ax[0].grid(alpha=.3); ax[0].axvline(0, color="k", lw=.5)
ax[0].axvline(-u10_eq, color="C1", ls=":", lw=1)   # equator value for reference
ax[1].plot(tau_x, lat, color="C3"); ax[1].set_xlabel(r"wind stress $\tau_x$ (N m$^{-2}$)   [$-$ = westward]")
ax[1].set_title(r"Zonal wind stress  ($C_D = 1.2\times10^{-3}$)"); ax[1].grid(alpha=.3); ax[1].axvline(0, color="k", lw=.5)
for a in ax: a.axhline(jet_lat, color="gray", ls=":"); a.axhline(-jet_lat, color="gray", ls=":")
save(fig, "02_wind_stress.png")
print(f"   τ_x at the equator = {rhoa*cD*u10_eq**2:.3f} N/m², peak = {np.max(np.abs(tau_x)):.3f} N/m² at ±{lat[np.argmax(np.abs(tau_x))]:.1f}°, wall = {rhoa*cD*speed[-1]**2:.3f} N/m² ({speed[-1]:.1f} m/s)")

# ================= 3. surface thermal forcing (pure SST restoring) =================
x = np.linspace(0, Lx, 600); xn = x/Lx
x_front, front_width = 0.6, 0.18
T_warm, T_cold, lam_T = 29.5, 24.0, 4e-5
warm = 0.5*(1-np.tanh((xn-x_front)/front_width))
T_target = T_cold + (T_warm-T_cold)*warm
fig, ax = plt.subplots(figsize=(8, 3.5))
ax.plot(x/1e3, T_target, color="C1"); ax.set_ylabel("°C"); ax.set_xlabel("x (km, W→E)")
ax.set_title(r"SST restoring target $T_{\rm target}(x)$,  $\lambda_T = 4\times10^{-5}$ m s$^{-1}$ (≈160 W m$^{-2}$ K$^{-1}$)")
ax.grid(alpha=.3)
save(fig, "03_surface_thermal_forcing.png")

# ================= 5. momentum sponge + deep drag =================
sw_x, sw_y = 5e5, 5e5
Xg = np.linspace(0, Lx, 200); Yg = np.linspace(-Ly/2, Ly/2, 150)
XX, YY = np.meshgrid(Xg, Yg)
west = np.clip((sw_x-XX)/sw_x, 0, None); east = np.clip((XX-(Lx-sw_x))/sw_x, 0, None)
southm = np.clip((sw_y-(YY+Ly/2))/sw_y, 0, None); northm = np.clip((YY-(Ly/2-sw_y))/sw_y, 0, None)
mask = np.maximum(west**2, east**2) + np.maximum(southm**2, northm**2)
zd = np.linspace(-1000, 0, 300)
z_top, z_bot = -300.0, -550.0
deepmask = np.clip((z_top - zd)/(z_top - z_bot), 0, 1)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5), gridspec_kw={"width_ratios":[1.6,1]})
pc = ax[0].pcolormesh(Xg/1e3, Yg/deg, mask, cmap="magma", shading="auto")
ax[0].set_xlabel("x (km)"); ax[0].set_ylabel("latitude (°)")
ax[0].set_title("Boundary sponge mask  (damps u, v)"); fig.colorbar(pc, ax=ax[0])
ax[1].plot(deepmask, zd, color="C4"); ax[1].set_xlabel("drag fraction")
ax[1].set_ylabel("z (m)"); ax[1].set_title("Deep Rayleigh drag  (confines EUC)")
ax[1].axhline(z_top, color="gray", ls=":"); ax[1].axhline(z_bot, color="gray", ls=":"); ax[1].grid(alpha=.3)
save(fig, "05_momentum_sponge_deep_drag.png")

# ================= 6. poleward heat sink =================
ys = np.linspace(-Ly/2, Ly/2, 600)
southb = np.clip((sw_y-(ys+Ly/2))/sw_y, 0, None); northb = np.clip((ys-(Ly/2-sw_y))/sw_y, 0, None)
nsband = np.maximum(southb**2, northb**2)
T_surf_ref, T_deep_ref, z_ref, d_ref = 22.0, 3.5, -150.0, 120.0
z_sink_mid, d_sink = -300.0, 75.0
zr = np.linspace(-1500, 0, 300)
Tref = T_deep_ref + (T_surf_ref-T_deep_ref)*0.5*(1+np.tanh((zr-z_ref)/d_ref))
zmask = 0.5*(1+np.tanh((zr-z_sink_mid)/d_sink))
# initial T profile for comparison
T_s0, T_mid, T_d0 = 28.0, 12.0, 3.5
zth, dth, zdt, ddt = -140.0, 55.0, -900.0, 450.0
Ti = T_d0 + (T_mid-T_d0)*0.5*(1+np.tanh((zr-zdt)/ddt)) + (T_s0-T_mid)*0.5*(1+np.tanh((zr-zth)/dth))
fig, ax = plt.subplots(1, 3, figsize=(14, 4.8))
ax[0].plot(ys/deg, nsband, color="C0"); ax[0].set_xlabel("latitude (°)")
ax[0].set_ylabel("sink strength"); ax[0].set_title("N/S sink band mask  $m_{NS}(y)$"); ax[0].grid(alpha=.3)
ax[1].plot(Tref, zr, label="extratropical reference $T_{ref}(z)$", color="C0")
ax[1].plot(Ti, zr, "--", label="initial profile $T_i(z)$", color="C1")
ax[1].set_xlabel("T (°C)"); ax[1].set_ylabel("z (m)"); ax[1].set_title("Heat-sink reference profile")
ax[1].legend(); ax[1].grid(alpha=.3)
ax[2].plot(zmask, zr, color="C2"); ax[2].set_xlabel("depth mask $m_z(z)$"); ax[2].set_ylabel("z (m)")
ax[2].set_title("Sink depth mask  ($z_m=-300$ m, $d_m=75$ m)"); ax[2].set_ylim(-1000, 0); ax[2].grid(alpha=.3)
save(fig, "06_poleward_heat_sink.png")

# ================= 7. initial conditions =================
S_s, S_sub, S_d = 34.6, 35.1, 34.8
zh, dh, zdf, ddf = -120.0, 70.0, -900.0, 500.0
Si = S_d + (S_sub-S_d)*0.5*(1+np.tanh((zr-zdf)/ddf)) + (S_s-S_sub)*0.5*(1+np.tanh((zr-zh)/dh))
fig, ax = plt.subplots(1, 2, figsize=(11, 5))
ax[0].plot(Ti, zr, color="C3"); ax[0].set_xlabel("T (°C)"); ax[0].set_ylabel("z (m)")
ax[0].set_title("Initial temperature $T_i(z)$"); ax[0].grid(alpha=.3)
ax[1].plot(Si, zr, color="C0"); ax[1].set_xlabel("S (g kg$^{-1}$)")
ax[1].set_title("Initial salinity $S_i(z)$"); ax[1].grid(alpha=.3)
save(fig, "07_initial_conditions.png")

# ================= 8. result snapshot (equatorial cross section) =================
# Candidate locations of the production spin-up's 3-D output (first one found wins).
RESULT_RUNS = [
    os.path.join(HERE, "..", "runs", "toy_model_spinup_2026-10-21_140641", "toy_model.jld2"),
    os.path.expanduser("~/Downloads/toy_model_spinup_2026-10-21_140641/toy_model.jld2"),
]
run = next((p for p in RESULT_RUNS if os.path.exists(p)), None)
if run:
    f = h5py.File(run, "r"); ts = f["timeseries"]; H = 4
    its = sorted(int(k) for k in ts["T"].keys() if k != "serialized")
    # last-100-day mean of the 30-day-mean snapshots (last 3 → days 1710–1800)
    sel = [str(k) for k in its[-3:]]
    days = [ts["t"][k][()]/86400 for k in sel]
    T = np.mean([np.array(ts["T"][k])[H:-H, H:-H, H:-H] for k in sel], axis=0)   # (z,y,x)
    u = np.mean([np.array(ts["u"][k])[H:-H, H:-H, H:-H] for k in sel], axis=0)
    yc = (np.arange(Ny)+0.5)*Ly/Ny - Ly/2
    band = np.abs(yc) <= 250e3
    Te = T[:, band, :].mean(axis=1); ue = u[:, band, :].mean(axis=1)
    if ue.shape[1] == Te.shape[1] + 1:        # u on x-faces → centres
        ue = 0.5*(ue[:, :-1] + ue[:, 1:])
    xc = (np.arange(Te.shape[1])+0.5)*Lx/Te.shape[1]
    fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    c0 = ax[0].contourf(xc/1e3, zc, Te, levels=20, cmap="turbo")
    ax[0].contour(xc/1e3, zc, Te, levels=[20.0], colors="w", linewidths=1.5)
    ax[0].set_ylim(-400, 0); ax[0].set_ylabel("z (m)"); fig.colorbar(c0, ax=ax[0], label="T (°C)")
    ax[0].set_title(f"Equatorial-band temperature (|y|≤250 km), days {days[0]-30:.0f}–{days[-1]:.0f} mean; white = 20 °C")
    m = np.nanpercentile(np.abs(ue), 99)
    c1 = ax[1].contourf(xc/1e3, zc, ue, levels=np.linspace(-m, m, 22), cmap="RdBu_r", extend="both")
    ax[1].set_ylim(-400, 0); ax[1].set_ylabel("z (m)"); ax[1].set_xlabel("x (km, W→E)")
    fig.colorbar(c1, ax=ax[1], label="u (m/s)")
    ax[1].set_title("Equatorial-band zonal velocity — eastward EUC (red) under westward SEC (blue)")
    save(fig, "08_result_equatorial_section.png")
    f.close()
else:
    print("spin-up output not found, skipping result snapshot (08)")

# ================= 9. warming experiments: forcing profiles =================
# Protocol A: ΔT_target added to both SST targets.  B: A plus ΔT_sink on the sink's surface
# reference.  C: warm pool +ΔT, cold tongue +2ΔT (ΔT_east = ΔT on T_cold only).
def target(dT, dT_east=0.0):
    return (T_cold + dT + dT_east) + ((T_warm + dT) - (T_cold + dT + dT_east))*warm
protocols = [("control",        0.0, 0.0, 0.0, "k",  "-"),
             ("A1  (+1 K)",     1.0, 0.0, 0.0, "C1", "-"),
             ("A2  (+2 K)",     2.0, 0.0, 0.0, "C3", "-"),
             ("A4  (+4 K)",     4.0, 0.0, 0.0, "C4", "-"),
             ("C2  (+2 W, +4 E)", 2.0, 2.0, 0.0, "C2", "--")]
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
for name, dT, dTe, dTs, col, ls in protocols:
    ax[0].plot(x/1e3, target(dT, dTe), color=col, ls=ls, label=name)
    ax[1].plot(x/1e3, target(dT, dTe) - target(0.0), color=col, ls=ls, label=name)
ax[1].plot(x/1e3, target(2.0) - target(0.0), color="C0", ls=":", lw=2.5, label="B2  (+2 K, surface = A2)")
ax[0].set_title("SST restoring target $T_{\\rm target}(x)$ by experiment"); ax[0].set_ylabel("°C")
ax[1].set_title("Imposed SST-target change, experiment − control"); ax[1].set_ylabel("ΔT (K)")
for a in ax[:2]: a.set_xlabel("x (km, W→E)"); a.grid(alpha=.3); a.legend(fontsize=8)
# sink reference: control vs B2
TrefB = (T_deep_ref) + ((T_surf_ref + 2.0)-T_deep_ref)*0.5*(1+np.tanh((zr-z_ref)/d_ref))
ax[2].plot(Tref, zr, color="k", label="control: $T_{\\rm surf}^{\\rm ref}=22$ °C")
ax[2].plot(TrefB, zr, color="C0", ls=":", lw=2.5, label="B2: $T_{\\rm surf}^{\\rm ref}=24$ °C")
ax[2].plot(Tref, zr, color="C3", ls="-", lw=0, label="A, C: unchanged")
ax[2].set_ylim(-800, 0); ax[2].set_xlabel("T (°C)"); ax[2].set_ylabel("z (m)")
ax[2].set_title("Sink reference profile $T_{\\rm ref}(z)$ (N/S bands)"); ax[2].grid(alpha=.3); ax[2].legend(fontsize=8)
save(fig, "09_warming_experiments.png")
print("done")
