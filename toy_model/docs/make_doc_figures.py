"""Generate boundary-condition / forcing figures for the toy-model description.
All curves are computed directly from the model's functional forms in toy_model.jl.
Figures are written to docs/figures/."""
import numpy as np
import matplotlib.pyplot as plt
import os, h5py

HERE = os.path.dirname(os.path.abspath(__file__))
FIG  = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

# ---- domain / constants (from toy_model.jl) ----
Lx, Ly, Lz = 1.5e7, 2.22e6, 2000.0
Nx, Ny, Nz = 140, 40, 60
rho0, cp, rho_fw, Lv = 1026.0, 3991.0, 1000.0, 2.5e6
deg = 111e3  # m per degree latitude

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
u10_trades, jet_lat, jet_width, cD, rhoa = 8.0, 5.0, 2.5, 2e-3, 1.2
y = np.linspace(-Ly/2, Ly/2, 400); lat = y/deg; lat_edge = (Ly/2)/deg
npk = np.exp(-((lat-jet_lat)/jet_width)**2)
spk = np.exp(-((lat+jet_lat)/jet_width)**2)
peak = np.maximum(npk, spk); eqs = np.exp(-(lat/8.0)**2)
shape = (peak + 0.9*eqs)*np.cos(np.pi*lat/(2*lat_edge))**2/1.2
u10 = u10_trades*shape
tau_kin = rhoa/rho0*cD*u10*np.abs(u10)     # m^2 s^-2
tau_Nm2 = tau_kin*rho0                       # N m^-2
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(u10, lat); ax[0].set_xlabel(r"$u_{10}$ (m s$^{-1}$)"); ax[0].set_ylabel("latitude (°)")
ax[0].set_title("10 m wind speed profile"); ax[0].grid(alpha=.3); ax[0].axhline(0, color="k", lw=.5)
ax[1].plot(tau_Nm2, lat, color="C3"); ax[1].set_xlabel(r"$\tau_x$ (N m$^{-2}$)")
ax[1].set_title("Zonal wind stress (trade winds)"); ax[1].grid(alpha=.3); ax[1].axhline(0, color="k", lw=.5)
for a in ax: a.axhline(jet_lat, color="gray", ls=":"); a.axhline(-jet_lat, color="gray", ls=":")
save(fig, "02_wind_stress.png")

# ================= 3. surface thermal forcing =================
x = np.linspace(0, Lx, 400); xn = x/Lx
x_front, front_width = 0.6, 0.18
T_warm, T_cold = 29.5, 24.0
warm = 0.5*(1-np.tanh((xn-x_front)/front_width))
T_target = T_cold + (T_warm-T_cold)*warm
E_base, E_warm = 1e-8, 3e-8
E = E_base + E_warm*warm
Q_latent = rho_fw*Lv*E       # W m^-2 (heat OUT)
Q_solar = 70.0
fig, ax = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
ax[0].plot(x/1e3, T_target, color="C1"); ax[0].set_ylabel("°C")
ax[0].set_title(r"SST restoring target $T_{\rm target}(x)$  (warm pool W → cold tongue E)")
ax[0].grid(alpha=.3)
ax[1].plot(x/1e3, E*1e3*86400, color="C2"); ax[1].set_ylabel("mm day$^{-1}$")
ax[1].set_title(r"Evaporation rate $E(x)$ (warm-pool enhanced)"); ax[1].grid(alpha=.3)
ax[2].plot(x/1e3, -Q_latent, color="seagreen", label="latent (evaporation) OUT")
ax[2].axhline(Q_solar, color="darkorange", label="solar IN (+70)")
ax[2].axhline(0, color="k", lw=.6)
ax[2].set_ylabel("W m$^{-2}$"); ax[2].set_xlabel("x (km, W→E)")
ax[2].set_title("Surface heat-flux components (into ocean +)"); ax[2].legend(); ax[2].grid(alpha=.3)
save(fig, "03_surface_thermal_forcing.png")

# (figure 04 "N/S boundary T flux" removed — those wall fluxes were deleted from
#  the model; the N/S walls are now insulating. Figure numbers 05–08 kept as-is.)

# ================= 5. momentum sponge + deep drag =================
sw_x, sw_y = 5e5, 5e5
Xg = np.linspace(0, Lx, 200); Yg = np.linspace(-Ly/2, Ly/2, 120)
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
ys = np.linspace(-Ly/2, Ly/2, 400)
southb = np.clip((sw_y-(ys+Ly/2))/sw_y, 0, None); northb = np.clip((ys-(Ly/2-sw_y))/sw_y, 0, None)
nsband = np.maximum(southb**2, northb**2)
T_surf_ref, T_deep_ref, z_ref, d_ref = 22.0, 3.5, -150.0, 120.0
zr = np.linspace(-1500, 0, 300)
Tref = T_deep_ref + (T_surf_ref-T_deep_ref)*0.5*(1+np.tanh((zr-z_ref)/d_ref))
# initial T profile for comparison
T_s0, T_mid, T_d0 = 28.0, 12.0, 3.5
zth, dth, zdt, ddt = -140.0, 55.0, -900.0, 450.0
Ti = T_d0 + (T_mid-T_d0)*0.5*(1+np.tanh((zr-zdt)/ddt)) + (T_s0-T_mid)*0.5*(1+np.tanh((zr-zth)/dth))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.8))
ax[0].plot(ys/deg, nsband, color="C0"); ax[0].set_xlabel("latitude (°)")
ax[0].set_ylabel("sink strength"); ax[0].set_title("N/S sink band mask  $m_{NS}(y)$"); ax[0].grid(alpha=.3)
ax[1].plot(Tref, zr, label="extratropical reference $T_{ref}(z)$", color="C0")
ax[1].plot(Ti, zr, "--", label="initial profile $T_i(z)$", color="C1")
ax[1].set_xlabel("T (°C)"); ax[1].set_ylabel("z (m)"); ax[1].set_title("Heat-sink reference profile")
ax[1].legend(); ax[1].grid(alpha=.3)
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
run = os.path.join(HERE, "..", "..", "runs", "toy_model_2026-06-05_234636", "toy_model.jld2")
if os.path.exists(run):
    f = h5py.File(run, "r"); ts = f["timeseries"]; H = 4
    its = sorted(int(k) for k in ts["t"].keys()); k = str(its[-1])
    day = ts["t"][k][()]/86400
    T = np.array(ts["T"][k])[H:-H, H:-H, H:-H]   # (z,y,x)
    u = np.array(ts["u"][k])[H:-H, H:-H, H:-H]
    yc = (np.arange(Ny)+0.5)*Ly/Ny - Ly/2
    band = np.abs(yc) <= 250e3
    Te = T[:, band, :].mean(axis=1); ue = u[:, band, :].mean(axis=1)
    if ue.shape[1] == Te.shape[1] + 1:        # u on x-faces → centres
        ue = 0.5*(ue[:, :-1] + ue[:, 1:])
    xc = (np.arange(Te.shape[1])+0.5)*Lx/Te.shape[1]
    fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    c0 = ax[0].contourf(xc/1e3, zc, Te, levels=20, cmap="turbo")
    ax[0].set_ylim(-400, 0); ax[0].set_ylabel("z (m)"); fig.colorbar(c0, ax=ax[0], label="T (°C)")
    ax[0].set_title(f"Equatorial-band temperature (|y|≤250 km), day {day:.0f}")
    m = np.nanpercentile(np.abs(ue), 99)
    c1 = ax[1].contourf(xc/1e3, zc, ue, levels=np.linspace(-m, m, 22), cmap="RdBu_r", extend="both")
    ax[1].set_ylim(-400, 0); ax[1].set_ylabel("z (m)"); ax[1].set_xlabel("x (km, W→E)")
    fig.colorbar(c1, ax=ax[1], label="u (m/s)")
    ax[1].set_title("Equatorial-band zonal velocity — eastward EUC (red) under westward SEC (blue)")
    save(fig, "08_result_equatorial_section.png")
else:
    print("run file not found, skipping result snapshot")
print("done")
