# Equatorial Pacific Toy Model — Model Description

A three-dimensional idealised ocean model of the equatorial Pacific, built with
[Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl). The goal is the
smallest configuration that reproduces the key features of the equatorial
Pacific:

1. a **tilted thermocline** (deep in the west, shallow in the east),
2. an **Equatorial Undercurrent (EUC)** — a subsurface eastward jet confined
   above ~300 m,
3. **little zonal momentum below the EUC**, and
4. **long-term stability** (a statistically steady state, not a drifting one).

It is driven by a constant trade-wind stress, a surface heat flux, a restoring
salinity flux, and a small set of body forcings that keep the solution physical
and bounded.

---

## 1. Domain and grid

The domain is a closed rectangular box on an equatorial β-plane.

| Quantity | Symbol | Value |
|---|---|---|
| Zonal extent | $L_x$ | $1.5\times10^7$ m (≈ 15 000 km) |
| Meridional extent | $L_y$ | $2.22\times10^6$ m (≈ ±10°) |
| Depth | $L_z$ | 25000 m |
| Grid points | $N_x\times N_y\times N_z$ | $140\times40\times60$ |
| Horizontal spacing | $\Delta x,\ \Delta y$ | ≈ 107 km, ≈ 55 km |
| Topology | — | bounded in all three directions |

The vertical grid is stretched: cells are uniform
(≈ 9.5 m) in the upper 300 m, and stretched till domain depth at 2000m. Cell thickness is given as:

$$
z_k = L_z\,\big(\zeta_0(k)\,\Sigma(k) - 1\big),\qquad
\zeta_0(k) = 1 + \frac{h(k)-1}{r},\quad
\Sigma(k) = \frac{1-e^{-s\,h(k)}}{1-e^{-s}},\quad h(k)=\tfrac{k-1}{N_z},
$$

with refinement $r=3.5$ and stretching $s=18$, modified from oceananigans example.

![Vertical grid spacing](figures/01_vertical_grid.png)

---

## 2. Governing equations

The model integrates the hydrostatic, and Boussinesqq, with a
free surface. For velocity $\mathbf{u}=(u,v,w)$, the horizontal momentum,
hydrostatic, and continuity equations are:

$$
\frac{\partial \mathbf{u}_h}{\partial t}
+ (\mathbf{u}\cdot\nabla)\,\mathbf{u}_h
+ f\,\hat{\mathbf{z}}\times\mathbf{u}_h
= -\nabla_h p + \nabla\cdot\boldsymbol{\tau} + \mathbf{F},
\qquad
\frac{\partial p}{\partial z} = b,
\qquad
\nabla\cdot\mathbf{u} = 0 .
$$

Conservative temperature $T$ and absolute salinity $S$ are given by:

$$
\frac{\partial T}{\partial t} + \mathbf{u}\cdot\nabla T = \nabla\cdot(\boldsymbol{\kappa}\nabla T) + F_T,
\qquad
\frac{\partial S}{\partial t} + \mathbf{u}\cdot\nabla S = \nabla\cdot(\boldsymbol{\kappa}\nabla S) .
$$

Equation of state is nonlinear TEOS-10, buoyancy is given as:

$$
b = -\,\frac{g\,\big(\rho(T,S,z) - \rho_0\big)}{\rho_0},\qquad \rho_0 = 1026~\mathrm{kg\,m^{-3}} .
$$

We use a $\beta$-plane approximation, $f = \beta\,y$, with
$\beta = 2\Omega/R = 2.29\times10^{-11}\ \mathrm{m^{-1}s^{-1}}$
($\Omega = 7.29\times10^{-5}\,\mathrm{s^{-1}}$, $R = 6.371\times10^{6}$ m). 

$\mathbf{F}$ and $F_T$ collect the body forcings (sponges, drag, heat sink)
described in §5.

---

## 3. Surface momentum forcing — trade winds

The wind is a constant easterly trade-wind stress with a realistic meridional
structure, ramped on over the first 10 days. A 10 m wind speed profile is built
from off-equatorial maxima at $\pm\varphi_j$ plus an equatorial enhancement:

$$
u_{10}(y) = \frac{u_{10}^{t}}{1.858}\,
\Big[\max\!\big(e^{-((\varphi-\varphi_j)/\Delta\varphi)^2},\,
e^{-((\varphi+\varphi_j)/\Delta\varphi)^2}\big)
+ 0.9\,e^{-(\varphi/8)^2}\Big],
$$

where $\varphi = y/(111\,\mathrm{km})$ is latitude in degrees. The profile is
**not** tapered to zero at the meridional walls: a zonal stress is *tangential*
to those walls (not a singularity), and the velocity sponge (§5.1) absorbs any
edge flow. Forcing the wind to zero there (an earlier `cos²` taper) imposed a
sharp $\partial\tau_x/\partial y$ — a spurious wind-stress curl — that drove
unphysical eastward subtropical flow near $\pm$10°. The normalization (1.858)
keeps the peak wind at $\approx$7 m s⁻¹; the residual at the wall is
$\approx$0.9 m s⁻¹, tangential and sponge-damped.

The stress is applied using Large & Yeager bulk formula with feed back from ocean surface currents. The atmospheric wind in an eastward-positive frame
($\mathbf{U}_a = (-u_{10},\,0)$ for easterly trades) and the ocean surface
velocity $\mathbf{u}_o = (u_o, v_o)$,

$$
\boldsymbol{\tau}(y,t) = \frac{\rho_a}{\rho_0}\,C_D\,
(\mathbf{U}_a - \mathbf{u}_o)\,\big|\mathbf{U}_a - \mathbf{u}_o\big|\;\cdot\;
\min\!\Big(1,\frac{t}{\tau_{\rm spin}}\Big),
\qquad
\big|\mathbf{U}_a - \mathbf{u}_s\big| = \sqrt{(u_{10}+u_o)^2 + v_o^2}.
$$

Both stress components are therefore boundary conditions that depend on the
surface velocity. 

| $u_{10}^{t}$ | $\varphi_j$ | $\Delta\varphi$ | $C_D$ | $\rho_a$ | $\tau_{\rm spin}$ |
|---|---|---|---|---|---|
| 8 m s⁻¹ | 5° | 2.5° | $2\times10^{-3}$ | 1.2 kg m⁻³ | 10 days |

![Wind stress profile](figures/02_wind_stress.png)

*The figure shows the resting-ocean reference stress (the relative-wind stress at
$\mathbf{u}_s=0$); the actual applied stress is reduced by the current feedback.*

(Sign convention: the stress enters as a flux of momentum through the top
boundary, and Oceananigans' top flux has the opposite sign to the physical
stress — a positive top $u$-flux removes eastward momentum.)


---

## 4. Surface tracer boundary conditions

### 4.1 Surface heat flux

The top temperature flux combines three physically distinct pieces. In the
Oceananigans sign convention a positive top flux removes heat:

$$
F_T^{\rm top} = \underbrace{\lambda_T\big(T - T_{\rm target}(x)\big)}_{\text{SST restoring}}
\;+\; \underbrace{\frac{\rho_{\rm fw}L_v}{\rho_0 c_p}\,E(x)}_{\text{latent (evaporative) cooling}}
\;-\; \underbrace{\frac{Q_{\rm solar}}{\rho_0 c_p}}_{\text{solar heating}} .
$$

**Restoring** is a *weak* relaxation toward a prescribed SST climatology with
piston velocity $\lambda_T = 7\times10^{-6}\,\mathrm{m\,s^{-1}}$
($\approx 29\ \mathrm{W\,m^{-2}\,K^{-1}}$, an air–sea damping timescale of ≈ 83 days
over a 50 m layer). The sign is a **negative feedback**: when $T>T_{\rm target}$
the flux is positive and cools. The target is zonally structured to set up the
warm-pool / cold-tongue contrast,

$$
T_{\rm target}(x) = T_{\rm cold} + (T_{\rm warm}-T_{\rm cold})\,
\tfrac12\Big[1 - \tanh\!\Big(\frac{x/L_x - x_f}{w_f}\Big)\Big].
$$

**Latent cooling** represents evaporation, enhanced in the warm pool, and is the
ocean's main surface heat *sink*. **Solar heating** is a uniform input that
roughly balances the area-mean evaporative cooling so the net surface budget is
≈ 0 and the SST stays near target.

| $\lambda_T$ | $T_{\rm warm}$ | $T_{\rm cold}$ | $x_f$ | $w_f$ | $Q_{\rm solar}$ |
|---|---|---|---|---|---|
| $7\times10^{-6}$ m s⁻¹ | 29.5 °C | 24.0 °C | 0.6 | 0.18 | 70 W m⁻² |

![Surface thermal forcing](figures/03_surface_thermal_forcing.png)

### 4.2 Evaporation and surface salinity flux

A single warm-pool-enhanced evaporation rate drives **both** the latent cooling
above and the surface salinity flux,

$$
E(x) = E_{\rm base} + E_{\rm warm}\,
\tfrac12\Big[1 - \tanh\!\Big(\frac{x/L_x - x_f}{w_f}\Big)\Big],
\qquad
F_S^{\rm top} = -E(x)\,S \;-\; r_S\,(S - S^{*}).
$$

The first term concentrates salt as freshwater leaves (evaporation makes the
surface saltier); the second is a weak relaxation of surface salinity toward a
target $S^{*}$.

| $E_{\rm base}$ | $E_{\rm warm}$ | $L_v$ | $\rho_{\rm fw}$ | $r_S$ | $S^{*}$ |
|---|---|---|---|---|---|
| $1\times10^{-8}$ m s⁻¹ | $3\times10^{-8}$ m s⁻¹ | $2.5\times10^{6}$ J kg⁻¹ | 1000 kg m⁻³ | $5\times10^{-7}$ m s⁻¹ | 34.6 g kg⁻¹ |

### 4.3 Meridional (north / south) walls

The north and south walls are **insulating** (zero heat flux), like the bottom
and the salinity walls. The meridional heat exchange and the cool subtropical
column are set by the poleward heat sink (§ 5.2).

### 4.4 Bottom

The bottom is insulating and free-slip.

---

## 5. Body forcings

### 5.1 Momentum sponges and deep Rayleigh drag

A linear, Rayleigh, damping is applied to the horizontal velocities,

$$
\mathbf{F}_{u,v} = -\big[\lambda_0\,m_{\rm sp}(x,y) + \lambda_d\,m_d(z)\big]\,(u,v).
$$

The boundary sponge $m_{\rm sp}$ damps velocities within 500 km of the
east/west and north/south walls (quadratic ramps), absorbing boundary reflections:

$$
m_{\rm sp}(x,y) = \alpha_{\rm EW}\max(w_W^2,w_E^2) + \alpha_{\rm NS}\max(w_S^2,w_N^2).
$$

The deep Rayleigh drag $m_d(z)$ is zero in the upper ocean (where the EUC
lives) and ramps to full strength below ≈ 550 m,

$$
m_d(z) = \mathrm{clamp}\!\Big(\frac{z_{\rm top}-z}{z_{\rm top}-z_{\rm bot}},\,0,\,1\Big),
\qquad z_{\rm top}=-300\,\mathrm{m},\ z_{\rm bot}=-550\,\mathrm{m}.
$$


While unphysical, this drag means that the EUC won't extend past 300m, and there is a removal of zonal monetum which would otherwise build up and mix everything.

| $\lambda_0$ | sponge width | $\lambda_d$ | $z_{\rm top}$ | $z_{\rm bot}$ |
|---|---|---|---|---|
| $1/5$ days | 500 km (each wall) | $1/20$ days | −300 m | −550 m |

![Momentum sponge and deep drag](figures/05_momentum_sponge_deep_drag.png)

### 5.2 Poleward heat sink (north/south temperature relaxation)

In a closed basin the surface heat taken up in the cold tongue (~+1.15 PW)
has nowhere to go and the ocean warms indefinitely. To mimic the real
tropical→extratropical heat export, $T$ is relaxed toward a **cool
extratropical reference profile** within the north/south bands:

$$
F_T^{\rm sink} = -\lambda_{\rm sink}\,m_{\rm NS}(y)\,\big(T - T_{\rm ref}(z)\big),
\qquad
T_{\rm ref}(z) = T_{\rm deep}^{\rm ref}
+ (T_{\rm surf}^{\rm ref}-T_{\rm deep}^{\rm ref})\,
\tfrac12\Big[1 + \tanh\!\Big(\frac{z - z_{\rm ref}}{d_{\rm ref}}\Big)\Big].
$$

Warm water carried to the subtropical edges is cooled there, exporting the
tropical heat uptake so the basin equilibrates. The deep reference matches the
initial deep temperature (3.5 °C) so the abyss is not nudged.

| $\lambda_{\rm sink}$ | $T_{\rm surf}^{\rm ref}$ | $T_{\rm deep}^{\rm ref}$ | $z_{\rm ref}$ | $d_{\rm ref}$ |
|---|---|---|---|---|
| $1/20$ days | 22 °C | 3.5 °C | −150 m | 120 m |

![Poleward heat sink](figures/06_poleward_heat_sink.png)

---

## 6. Subgrid closure and advection

Vertical mixing combines a tiny background diffusivity (to preserve the
sharp thermocline) with convective adjustment — a large diffusivity that
fires only where the column is statically unstable ($N^2<0$). A hydrostatic
model cannot represent convection, so this acts as a safety valve that removes
static instabilities instead of letting them blow up:

$$
\kappa_z = \begin{cases}\kappa_{\rm conv} & N^2 < 0\\[2pt]\kappa_{\rm bg} & N^2 \ge 0\end{cases},
\qquad \kappa_{\rm conv}=1.0,\ \kappa_{\rm bg}=10^{-5},\ \nu_{\rm bg}=10^{-3}\ \mathrm{m^2\,s^{-1}} .
$$

**Horizontal mixing** is a **biharmonic viscosity on momentum only**,
$\partial_t\mathbf{u}_h \supset -\nu_4\,\nabla_h^4\mathbf{u}_h$ with
$\nu_4 = 10^{12}\ \mathrm{m^4\,s^{-1}}$. Being scale-selective, it damps the
grid-scale ($2\Delta x$) noise that would otherwise appear as vertical stripes
in $w$, while leaving the basin-scale EUC essentially untouched. It is *not*
applied to tracers, because biharmonic diffusion is not sign-preserving and can
drive salinity negative.

The vertical mixing uses an **implicit** time discretisation, so the large
convective diffusivity ($\kappa_{\rm conv}=1\ \mathrm{m^2\,s^{-1}}$) does not
constrain the time step. Convective adjustment is applied to tracers only
($\nu_{\rm conv}=0$), leaving the EUC momentum to the biharmonic viscosity and
the Rayleigh drag.

**Advection** uses WENO schemes: order 3 for momentum, and **order 7 for
tracers**. The high-order tracer scheme has far less spurious diapycnal mixing,
which keeps the thermocline (and the eastern cold-tongue source water) sharp.

Together these three choices are what make long integrations stable and
physical: the convective adjustment is the safety valve a hydrostatic model
*needs* (with no vertical momentum equation it cannot represent convection, so
without this term static instabilities grow without bound); the momentum-only
biharmonic removes the grid-scale velocity noise that would otherwise be
amplified — via spurious vertical motions — into column-wide mixing; and the
high-order tracer advection stops the thermocline from being eroded
numerically. Each was added in response to a specific failure mode (runaway
convection, checkerboard noise, and thermocline erosion, respectively).

---

## 7. Initial conditions

Both tracers start horizontally uniform with idealised tropical stratification:
a warm mixed layer over a sharp thermocline (temperature), and a subsurface
salinity maximum (halocline). Velocities start from rest.

$$
T_i(z) = T_d + (T_m-T_d)\,\tfrac12\!\left[1+\tanh\tfrac{z-z_{dt}}{\delta_{dt}}\right]
+ (T_s-T_m)\,\tfrac12\!\left[1+\tanh\tfrac{z-z_{th}}{\delta_{th}}\right]
$$

with $T_s=28$, $T_m=12$, $T_d=3.5$ °C, thermocline at $z_{th}=-140$ m
($\delta_{th}=55$ m). The salinity profile $S_i(z)$ has an analogous two-tanh
form ($S=34.6$ surface, 35.1 subsurface, 34.8 deep).

![Initial conditions](figures/07_initial_conditions.png)

---

## 8. Time stepping, output and diagnostics

The simulation uses an adaptive time step (CFL = 0.2, capped at 600 s). Full
3-D fields ($T,S,u,v,w,\partial_z T$) are saved every 48 h; scalar diagnostics
every 24 h. Two diagnostic families are computed online:

- **Momentum spin-down**: volume-integrated zonal momentum split into the
  upper ocean ($z\ge-300$ m) and the deep ocean ($z<-300$ m), including
  $\int|u|\,\mathrm{d}V$ and $\int u^2\,\mathrm{d}V$ below 300 m. The deep
  integrals flattening near zero confirms the deep ocean is spinning down.
- **Surface heat budget** (post-processed): the surface flux is reconstructed
  by mechanism (solar in, restoring in/out, evaporation out) and compared with
  the actual ocean heat-content tendency $\mathrm{d}/\mathrm{d}t\int\rho c_p T\,\mathrm{d}V$.
  The tendency approaching zero confirms thermal equilibration; the gap between
  it and the surface net flux is the poleward export through the heat sink.

---

## 9. Spin-up and equilibration

The model starts from rest with horizontally-uniform stratification (§7), so the
circulation has to spin up. The trade-wind stress is ramped on linearly over the
first 10 days,
$\tau_x(t)\propto\min(1,\,t/\tau_{\rm spin})$ with $\tau_{\rm spin}=10$ days,
and is steady thereafter. Two very different timescales then govern the approach
to a statistically steady state:

**Dynamical adjustment — fast (weeks to a few months).** Once the winds are up,
the surface flow and an equatorial wave field (Kelvin and Rossby modes)
redistribute mass: warm water piles up in the west, the thermocline tilts, and
the resulting subsurface zonal pressure gradient drives the EUC. A
first-baroclinic equatorial Kelvin wave crosses the basin in $\sim$70 days, so
the thermocline tilt and the EUC reach a recognisable, quasi-steady structure
within roughly the first 1–3 months.

**Thermal adjustment — slow (many months).** The *heat content* equilibrates
much more slowly, because the surface fluxes (solar, evaporation, restoring) and
the poleward heat sink must come into balance through the whole water column.
This is monitored with the heat-content tendency
$\mathrm{d}/\mathrm{d}t\!\int\!\rho c_p T\,\mathrm{d}V$ (§8): the run is thermally
equilibrated when this tendency settles to $\approx 0$. Likewise, the
deep-momentum diagnostic flattening near zero confirms the deep ocean has spun
down rather than accumulating zonal momentum.

Because of the slow thermal branch, the model is integrated for **300+ days**,
and time-mean diagnostics are formed over the **last 100 days** so they
represent the equilibrated state rather than the spin-up transient.

---

## 10. Resulting state

With this configuration the model reproduces the target structure: a thermocline
that tilts from ~225 m deep in the west to ~30 m in the east, an eastward EUC
core of ~0.5–0.7 m s⁻¹ near the thermocline beneath a westward surface current
(South Equatorial Current), and quiescent water below ~300 m.

![Equatorial cross section](figures/08_result_equatorial_section.png)

*Equatorial-band ($|y|\le250$ km) mean from a 300-day run. Top: temperature,
showing the tilted thermocline (warm west, cold-tongue outcrop in the east).
Bottom: zonal velocity, showing the eastward Equatorial Undercurrent (red,
core ≈ 100–150 m) beneath the westward surface flow (blue), and weak flow below.*

---

### Parameter summary

| Group | Parameters |
|---|---|
| Grid | $L_x=15{,}000$ km, $L_y=\pm10°$, $L_z=2000$ m; $140\times40\times60$ |
| Rotation | β-plane, $\beta=2.29\times10^{-11}$ m⁻¹s⁻¹ |
| Wind | $u_{10}^t=8$ m s⁻¹, peaks at ±5°, $C_D=2\times10^{-3}$, 10-day ramp |
| Heat flux | $\lambda_T=7\times10^{-6}$ m s⁻¹; $T_{\rm target}$: 29.5→24 °C; $Q_{\rm solar}=70$ W m⁻² |
| Evaporation | $E=(1\text{–}4)\times10^{-8}$ m s⁻¹ (warm-pool enhanced) |
| Salinity | restore $5\times10^{-7}$ m s⁻¹ to 34.6 g kg⁻¹ |
| Momentum drag | sponge $1/5$ day (walls); deep drag $1/20$ day below 300–550 m |
| Heat sink | $1/20$ day toward 22 °C in N/S bands |
| Closure | convective adjustment ($\kappa_{\rm conv}=1$) + biharmonic $\nu_4=10^{12}$ m⁴s⁻¹ |
| Advection | WENO-3 momentum, WENO-7 tracers |
