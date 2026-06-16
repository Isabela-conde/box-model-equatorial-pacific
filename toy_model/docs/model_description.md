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

It is driven by a constant trade-wind stress, a surface heat flux (restoring +
evaporative cooling), a restoring salinity flux, penetrating solar heating, and a
small set of body forcings that keep the solution physical and bounded.

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

$\mathbf{F}$ and $F_T$ collect the body forcings (momentum sponges, deep drag,
penetrating solar heating, optional poleward heat sink) described in §5.

---

## 3. Surface momentum forcing — trade winds

The wind is a constant easterly trade-wind stress with a realistic meridional
structure, ramped on over the first 10 days. The 10 m easterly wind speed is a
weak equatorial value $u_{10}^{\rm eq}$ (the doldrums minimum) plus a smooth
off-equatorial trade bump that peaks at $\pm\varphi_j$, minus an edge taper that
weakens the easterlies toward the north/south walls:

$$
u_{10}(y) = u_{10}^{\rm eq}
+ \big(u_{10}^{t} - u_{10}^{\rm eq}\big)\,
\underbrace{\Big(\tfrac{\varphi}{\varphi_j}\Big)^2 e^{\,1-(\varphi/\varphi_j)^2}}_{\text{trade bump}}
\;-\;
u_{10}^{\rm wall}\,
\underbrace{\Big(\tfrac{\varphi}{\varphi_e}\Big)^4}_{\text{wall taper}},
$$

where $\varphi = y/(111\,\mathrm{km})$ is latitude in degrees and
$\varphi_e = L_y/2/(111\,\mathrm{km}) \approx 10°$ is the wall latitude. The trade
bump is an even, $C^\infty$ function that is exactly $0$ at the equator (so
$u_{10} = u_{10}^{\rm eq}$ there) and exactly $1$ at $\pm\varphi_j$ (so the peak
is $u_{10}^{t}$); it replaces an earlier `max`-of-two-Gaussians form that had a
non-differentiable **cusp** at the equator. The wall taper $(\varphi/\varphi_e)^4$
is $\approx0$ through the interior — leaving the equatorial dip and the trade peak
essentially untouched — and reaches $1$ at the walls, where $u_{10}^{\rm wall}$ is
removed. The net shape: weakest-but-nonzero easterlies dipping at the equator
($u_{10}^{\rm eq}=5$ m s⁻¹), a trade maximum ($\approx5.9$ m s⁻¹) near
$\pm\varphi_j=\pm6°$, and the easterlies tapering to $\approx4$ m s⁻¹ at the walls
— **below** the equatorial value. The wind is deliberately **not** forced to zero
at the walls: a zonal stress is *tangential* there (not a singularity), the
velocity sponge (§5.1) absorbs edge flow, and a hard taper-to-zero would impose a
spurious wind-stress curl that drives unphysical eastward subtropical flow.

| $u_{10}^{\rm eq}$ | $u_{10}^{t}$ | $\varphi_j$ | $u_{10}^{\rm wall}$ |
|---|---|---|---|
| 5 m s⁻¹ | 6 m s⁻¹ | 6° | 1.5 m s⁻¹ (drop at walls) |

![Wind stress profile](figures/02_wind_stress.png)

*Surface zonal wind (left) and the resulting wind stress (right). The dotted line
marks the equatorial value; note the trade peak near $\pm6°$ and the walls
tapering below the equator.*

**Resting-ocean vs. relative-wind stress (a build-time switch).** The stress is a
bulk formula; the model can apply it two ways, selected by the `relative_wind`
keyword to `build_model` (and exposed through `run_model`):

- `relative_wind = false` (**resting-ocean / "ocean at rest", the default for the
  current runs**): the stress uses the wind speed alone,
  $\tau_x = \tfrac{\rho_a}{\rho_0} C_D\,u_{10}\,|u_{10}|$, with no top stress on
  $v$.
- `relative_wind = true` (**current feedback**): the stress uses the wind
  *relative* to the ocean surface current, with the atmospheric wind in an
  eastward-positive frame ($\mathbf{U}_a = (-u_{10},\,0)$ for easterly trades) and
  the ocean surface velocity $\mathbf{u}_o = (u_o, v_o)$,
  $$
  \boldsymbol{\tau}(y,t) = \frac{\rho_a}{\rho_0}\,C_D\,
  (\mathbf{U}_a - \mathbf{u}_o)\,\big|\mathbf{U}_a - \mathbf{u}_o\big|\;\cdot\;
  \min\!\Big(1,\frac{t}{\tau_{\rm spin}}\Big),
  $$
  so both stress components become boundary conditions that depend on the surface
  velocity (and $\tau_y$ opposes the ocean's own $v$, damping eddies). At rest the
  two forms are identical.

Both forms are ramped on linearly over $\tau_{\rm spin}=10$ days.

| $C_D$ | $\rho_a$ | $\tau_{\rm spin}$ |
|---|---|---|
| $2\times10^{-3}$ | 1.2 kg m⁻³ | 10 days |

(Sign convention: the stress enters as a flux of momentum through the top
boundary, and Oceananigans' top flux has the opposite sign to the physical
stress — a positive top $u$-flux removes eastward momentum.)


---

## 4. Surface tracer boundary conditions

### 4.1 Surface heat flux

The top temperature flux combines two surface pieces. In the Oceananigans sign
convention a positive top flux removes heat:

$$
F_T^{\rm top} = \underbrace{\lambda_T\big(T - T_{\rm target}(x)\big)}_{\text{SST restoring}}
\;+\; \underbrace{\frac{\rho_{\rm fw}L_v}{\rho_0 c_p}\,E(x)}_{\text{latent (evaporative) cooling}} .
$$

Solar heating is **not** a surface flux — it is applied as a depth-penetrating
*body* forcing (§5.3), because a surface flux deposits all the shortwave heat in
the single top grid cell, building a razor-sharp near-surface stratification that
caps the wind-mixed layer at $\sim$15 m and strangles the wind-driven SEC.

**Restoring** is a relaxation toward a prescribed SST climatology with
piston velocity $\lambda_T = 2\times10^{-5}\,\mathrm{m\,s^{-1}}$
($\approx 82\ \mathrm{W\,m^{-2}\,K^{-1}}$, an air–sea damping timescale of ≈ 29 days
over a 50 m layer). The sign is a **negative feedback**: when $T>T_{\rm target}$
the flux is positive and cools. The target is zonally structured to set up the
warm-pool / cold-tongue contrast,

$$
T_{\rm target}(x) = T_{\rm cold} + (T_{\rm warm}-T_{\rm cold})\,
\tfrac12\Big[1 - \tanh\!\Big(\frac{x/L_x - x_f}{w_f}\Big)\Big].
$$

**Latent cooling** represents evaporation, enhanced in the warm pool, and is the
ocean's main surface heat *sink*. It is balanced (in the area mean) by the
penetrating solar input (§5.3) so the net column heat budget is ≈ 0.

| $\lambda_T$ | $T_{\rm warm}$ | $T_{\rm cold}$ | $x_f$ | $w_f$ |
|---|---|---|---|---|
| $2\times10^{-5}$ m s⁻¹ | 29.5 °C | 24.0 °C | 0.6 | 0.18 |

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

The sink is a **build-time switch**: the `poleward_sink` keyword to `build_model`
(threaded through `run_model`) multiplies $F_T^{\rm sink}$ by 1 when on and 0 when
off. It is currently **off** — with the sink active the basin equilibrates but the
N/S bands cool the whole column; with it off the basin warms slowly but the
equatorial structure is cleaner. (The solar penetration of §5.3 is independent of
this switch and always on.)

| $\lambda_{\rm sink}$ | $T_{\rm surf}^{\rm ref}$ | $T_{\rm deep}^{\rm ref}$ | $z_{\rm ref}$ | $d_{\rm ref}$ |
|---|---|---|---|---|
| $1/150$ days | 22 °C | 3.5 °C | −150 m | 120 m |

![Poleward heat sink](figures/06_poleward_heat_sink.png)

### 5.3 Penetrating solar heating

Shortwave heating enters as a body forcing with an exponential penetration profile
rather than a surface flux. The downward solar flux decays with depth as
$I(z) = Q_{\rm solar}\,e^{z/h_{\rm solar}}$ ($z\le0$), and the local heating rate
is its convergence:

$$
F_T^{\rm solar}(z) = \frac{1}{\rho_0 c_p}\,\frac{\partial I}{\partial z}
= \frac{Q_{\rm solar}}{\rho_0 c_p\,h_{\rm solar}}\,e^{\,z/h_{\rm solar}}.
$$

Integrated over the column this absorbs exactly $Q_{\rm solar}$ — the same total
heat as the old surface term, so the area-mean budget still balances the latent
cooling — but spread over an $e$-folding depth $h_{\rm solar}$ instead of the top
cell. The top $\approx$9.5 m cell now receives only $1-e^{-\Delta z/h_{\rm solar}}
\approx 32\%$ of the shortwave, and the rest warms deeper water, so the wind-mixed
layer stays deep enough for the trade winds to drive a proper westward SEC.
Larger $h_{\rm solar}$ → deeper, gentler heating → deeper mixed layer → stronger
SEC.

| $Q_{\rm solar}$ | $h_{\rm solar}$ |
|---|---|
| 70 W m⁻² | 25 m |

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
- **Column heat budget** (post-processed): the heating is reconstructed by
  mechanism (solar in, restoring in/out, evaporation out — note solar now enters
  as a penetrating body forcing, §5.3, but its column integral is unchanged) and
  compared with the actual ocean heat-content tendency
  $\mathrm{d}/\mathrm{d}t\int\rho c_p T\,\mathrm{d}V$. The tendency approaching
  zero confirms thermal equilibration; when the poleward heat sink (§5.2) is
  enabled, the gap between the tendency and the surface net flux is the poleward
  export. With the sink off the basin instead drifts slowly warm.

---

## 9. Spin-up and equilibration

The model starts from rest with horizontally-uniform stratification (§7), so the
circulation has to spin up. The trade-wind stress is ramped on linearly over the
first 10 days,
$\tau_x(t)\propto\min(1,\,t/\tau_{\rm spin})$ with $\tau_{\rm spin}=10$ days,
and is steady thereafter. Two very different timescales then govern the approach
to a statistically steady state:

**Dynamical adjustment — fast (weeks to a few months), but with a large
overshoot.** Once the winds are up, the surface flow and an equatorial wave field
(Kelvin and Rossby modes) redistribute mass: warm water piles up in the west, the
thermocline tilts, and the resulting subsurface zonal pressure gradient drives the
EUC. A first-baroclinic equatorial Kelvin wave crosses the basin in $\sim$70 days.

Crucially, this adjustment is **not monotonic** — it overshoots. The typical
equatorial-band trajectory is:

- **Days ~10–60:** the winds drive a strong, correct westward surface SEC
  ($u_{\rm surf}\approx-1.5\,\mathrm{m\,s^{-1}}$) with a weak, deep EUC.
- **Days ~100–200:** the zonal pressure gradient builds faster than the basin can
  equilibrate, so the EUC transiently **shoals to $\sim$25 m and balloons**
  ($>1\,\mathrm{m\,s^{-1}}$), filling the upper layer with eastward flow and
  **nearly flipping the surface eastward** — the run *looks broken* here.
- **Days ~250–350:** the wave field reflects and settles; the EUC sinks back to
  the thermocline and the westward SEC re-establishes on top.

This overshoot is a genuine feature of the spin-up, **not** a configuration bug:
the reference run reproduces exactly this trajectory and only looks correct after
$\sim$day 250. **A run stopped at 150–200 days is caught mid-overshoot and will
show a shallow eastward "slab" with little or no SEC** — the fix is simply to
integrate longer, not to change the forcing.

**Thermal adjustment — slow (many months).** The *heat content* equilibrates
much more slowly, because the surface fluxes (solar, evaporation, restoring) and
the poleward heat sink must come into balance through the whole water column.
This is monitored with the heat-content tendency
$\mathrm{d}/\mathrm{d}t\!\int\!\rho c_p T\,\mathrm{d}V$ (§8): the run is thermally
equilibrated when this tendency settles to $\approx 0$. Likewise, the
deep-momentum diagnostic flattening near zero confirms the deep ocean has spun
down rather than accumulating zonal momentum.

Because of both the dynamical overshoot and the slow thermal branch, the model is
integrated for **≥300–400 days**, and time-mean diagnostics are formed over the
**last 100 days** so they represent the equilibrated state rather than the
spin-up transient. Judging the EUC/SEC from an earlier snapshot is misleading.

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
| Wind | $u_{10}^{\rm eq}=5$, $u_{10}^t=6$ m s⁻¹ (peak ±6°), wall drop 1.5 m s⁻¹; $C_D=2\times10^{-3}$, 10-day ramp; resting-ocean stress (`relative_wind=false`) |
| Heat flux | $\lambda_T=2\times10^{-5}$ m s⁻¹; $T_{\rm target}$: 29.5→24 °C (surface restoring + latent only) |
| Solar | $Q_{\rm solar}=70$ W m⁻², penetrating body forcing, $h_{\rm solar}=25$ m |
| Evaporation | $E=(1\text{–}4)\times10^{-8}$ m s⁻¹ (warm-pool enhanced) |
| Salinity | restore $5\times10^{-7}$ m s⁻¹ to 34.6 g kg⁻¹ |
| Momentum drag | sponge $1/5$ day (walls); deep drag $1/20$ day below 300–550 m |
| Heat sink | $1/150$ day toward 22 °C in N/S bands (`poleward_sink`, currently **off**) |
| Closure | convective adjustment ($\kappa_{\rm conv}=1$) + biharmonic $\nu_4=10^{12}$ m⁴s⁻¹ |
| Advection | WENO-3 momentum, WENO-7 tracers |
