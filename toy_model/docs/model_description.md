# Equatorial Pacific Toy Model — Model Description

A three-dimensional idealised ocean model of the equatorial Pacific, built with
[Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) and run on a single
GPU (NCI Gadi, V100). The goal is the smallest configuration that reproduces the
key features of the equatorial Pacific:

1. a **tilted thermocline** (deep in the west, shallow in the east),
2. an **Equatorial Undercurrent (EUC)** — a subsurface eastward jet confined
   above ~300 m, with a westward South Equatorial Current (SEC) on top,
3. **little zonal momentum below the EUC**, and
4. **long-term stability** (a statistically steady state, not a drifting one),

and then to use it for **surface-warming experiments** (Paper 1): how the EUC's
strength and depth respond when the sea surface is warmed (§11).

It is driven by a constant trade-wind stress, a pure SST-restoring surface heat
flux, an evaporation + restoring surface salinity flux, and a small set of body
forcings (momentum sponges, deep drag, a poleward heat sink) that keep the solution
physical and bounded. The configuration below is the **production configuration of
2026-10-21**; everything is computed from the functional forms in `toy_model.jl`.

---

## 1. Domain and grid

The domain is a closed rectangular box on an equatorial β-plane.

| Quantity | Symbol | Value |
|---|---|---|
| Zonal extent | $L_x$ | $1.5\times10^7$ m (≈ 15 000 km, read as 120°E → 80°W) |
| Meridional extent | $L_y$ | $3.33\times10^6$ m (≈ ±15°) |
| Depth | $L_z$ | 2000 m |
| Grid points | $N_x\times N_y\times N_z$ | $600\times150\times60$ |
| Horizontal spacing | $\Delta x,\ \Delta y$ | ≈ 25 km, ≈ 22 km |
| Halo | — | 4 cells (WENO) |
| Topology | — | bounded in all three directions |
| Architecture | — | `GPU()` |

The horizontal resolution (25 km) resolves the equatorial deformation radius and
tropical instability waves; at the earlier 107 km the eastern cold-tongue front
developed a standing $2\Delta x$ checkerboard that the closures could not remove.

The vertical grid is stretched: cells are uniform (≈ 9.5 m) in the upper 300 m and
stretch to ≈ 370 m at the 2000 m bottom:

$$
z_k = L_z\,\big(\zeta_0(k)\,\Sigma(k) - 1\big),\qquad
\zeta_0(k) = 1 + \frac{h(k)-1}{r},\quad
\Sigma(k) = \frac{1-e^{-s\,h(k)}}{1-e^{-s}},\quad h(k)=\tfrac{k-1}{N_z},
$$

with refinement $r=3.5$ and stretching $s=18$ (modified from an Oceananigans example).

![Vertical grid spacing](figures/01_vertical_grid.png)

---

## 2. Governing equations

The model integrates the hydrostatic Boussinesq equations with an implicit free
surface. For velocity $\mathbf{u}=(u,v,w)$, the horizontal momentum, hydrostatic,
and continuity equations are:

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

Conservative temperature $T$ and absolute salinity $S$ obey

$$
\frac{\partial T}{\partial t} + \mathbf{u}\cdot\nabla T = \nabla\cdot(\boldsymbol{\kappa}\nabla T) + F_T,
\qquad
\frac{\partial S}{\partial t} + \mathbf{u}\cdot\nabla S = \nabla\cdot(\boldsymbol{\kappa}\nabla S) .
$$

The equation of state is the nonlinear TEOS-10 polynomial; buoyancy is

$$
b = -\,\frac{g\,\big(\rho(T,S,z) - \rho_0\big)}{\rho_0},\qquad \rho_0 = 1026~\mathrm{kg\,m^{-3}} .
$$

Rotation is a $\beta$-plane, $f = \beta\,y$, with
$\beta = 2\Omega/R = 2.29\times10^{-11}\ \mathrm{m^{-1}s^{-1}}$
($\Omega = 7.29\times10^{-5}\,\mathrm{s^{-1}}$, $R = 6.371\times10^{6}$ m).

$\mathbf{F}$ and $F_T$ collect the body forcings (momentum sponges, deep drag,
poleward heat sink) described in §5.

---

## 3. Surface momentum forcing — trade winds

The wind is a constant easterly trade-wind stress with a smooth meridional
structure, ramped on over the first 10 days. The 10 m easterly wind speed has an
equatorial dip, off-equatorial trade maxima, and a quartic wall taper that weakens
(but does not remove) the easterlies at the north/south walls:

$$
u_{10}(\varphi) = u_{10}^{\rm eq}
+ \big(u_{10}^{\rm tr} - u_{10}^{\rm eq}\big)\, r^2 e^{\,1-r^2}
- u_{10}^{\rm drop}\Big(\frac{\varphi}{\varphi_e}\Big)^{4},
\qquad r = \frac{\varphi}{\varphi_j},
$$

where $\varphi = y/(111\,\mathrm{km})$ is latitude in degrees and
$\varphi_e = 15°$ is the wall latitude. The bump $r^2e^{1-r^2}$ is exactly 0 at the
equator and 1 at $\pm\varphi_j$, with no kink. The net shape: 7 m s⁻¹ at the
equator, a trade maximum of ≈ 7.7 m s⁻¹ near $\pm7°$ (the quartic term pulls the
peak equatorward of $\varphi_j$), and ≈ 3.6 m s⁻¹ at the walls.

| $u_{10}^{\rm eq}$ | $u_{10}^{\rm tr}$ | $\varphi_j$ | $u_{10}^{\rm drop}$ | equator | peak | wall |
|---|---|---|---|---|---|---|
| 7 m s⁻¹ | 8 m s⁻¹ | 10° | 4 m s⁻¹ | 7.0 m s⁻¹ | ≈ 7.7 m s⁻¹ at ±7° | ≈ 3.6 m s⁻¹ |

The stress is a bulk formula on the resting ocean (`relative_wind = false`):

$$
\tau_x = \frac{\rho_a}{\rho_0}\,C_D\,u_{10}\,|u_{10}|\;\min\!\Big(1,\frac{t}{\tau_{\rm spin}}\Big),
$$

| $C_D$ | $\rho_a$ | $\tau_{\rm spin}$ | $\tau_x$ equator | $\tau_x$ peak |
|---|---|---|---|---|
| $1.2\times10^{-3}$ | 1.2 kg m⁻³ | 10 days | 0.071 N m⁻² | ≈ 0.085 N m⁻² |

![Wind stress profile](figures/02_wind_stress.png)

*Surface zonal wind (left) and the resulting wind stress (right).*

> **Why $C_D = 1.2\times10^{-3}$.** The earlier value $2\times10^{-3}$ gave an
> equatorial stress of 0.118 N m⁻², about twice the observed annual mean at 140°W.
> In one-year branch experiments from a spun-up state it was the **only** change
> among several (tracer diffusivity, viscosity, drag) that restored a westward
> surface layer over the EUC, and it also reduced the cold-tongue bias from 6 K to
> about 2.5 K during the branch. The standard bulk value was therefore adopted and
> the model re-spun-up.

A current-feedback (relative-wind) form of the stress exists as a build-time switch
(`relative_wind = true`) but is **off** in all production runs.

(Sign convention: the stress enters as a flux of momentum through the top boundary,
and Oceananigans' top flux has the opposite sign to the physical stress — a positive
top $u$-flux removes eastward momentum.)

---

## 4. Surface tracer boundary conditions

### 4.1 Surface heat flux

The top temperature flux is **pure SST restoring**. In the Oceananigans sign
convention a positive top flux removes heat:

$$
F_T^{\rm top} = \lambda_T\big(T - T_{\rm target}(x)\big).
$$

It is a relaxation toward a prescribed SST climatology with piston velocity
$\lambda_T = 4\times10^{-5}\,\mathrm{m\,s^{-1}}$
($\approx 164\ \mathrm{W\,m^{-2}\,K^{-1}}$; an e-folding of ≈ 14 days over a 50 m
layer). The sign is a **negative feedback**: when $T>T_{\rm target}$ the flux is
positive and cools. The target is zonally structured to set up the warm-pool /
cold-tongue contrast,

$$
T_{\rm target}(x) = T_{\rm cold} + (T_{\rm warm}-T_{\rm cold})\,
\tfrac12\Big[1 - \tanh\!\Big(\frac{x/L_x - x_f}{w_f}\Big)\Big].
$$

| $\lambda_T$ | $T_{\rm warm}$ | $T_{\rm cold}$ | $x_f$ | $w_f$ |
|---|---|---|---|---|
| $4\times10^{-5}$ m s⁻¹ | 29.5 °C | 24.0 °C | 0.6 | 0.18 |

![Surface thermal forcing](figures/03_surface_thermal_forcing.png)

> **On the restoring strength.** $\lambda_T$ was doubled from $2\times10^{-5}$ on
> 2026-10-20 because the eastern cold tongue sat ~8 K below its target: the heat
> the restoring pumps into the cold tongue (~4 PW) is set by the equatorial
> upwelling, and $\lambda_T$ only sets how far below target the SST sits to deliver
> it. Doubling halved the deficit. The price is a stiffer SST, which damps the SST
> response in the warming experiments; warming is therefore imposed through the
> target (§11), not through a flux.

### 4.2 Surface salinity flux

A spatially **uniform** evaporation plus a weak relaxation toward a target surface
salinity (positive top flux removes salt):

$$
F_S^{\rm top} = -E\,S \;+\; r_S\,(S - S^{*}).
$$

The first term concentrates salt as freshwater leaves; the second relaxes surface
salinity toward $S^{*}$. The two balance at a surface salinity of
$S = r_S S^*/(r_S - E) \approx 35.3$, close to the equilibrated value (35.0).

| $E$ | $r_S$ | $S^{*}$ |
|---|---|---|
| $1\times10^{-8}$ m s⁻¹ | $5\times10^{-7}$ m s⁻¹ | 34.6 g kg⁻¹ |

> **History note.** Until 2026-10-21 the restoring term carried a minus sign
> (anti-restoring): surface salinity drifted 34.6 → 35.3 over a five-year run and
> never approached the target. The sign is now the same as the temperature
> restoring. The drift contributed only a few per cent to the upper-ocean
> stratification, so earlier dynamical conclusions stand, but every production run
> uses the corrected form.

### 4.3 Meridional walls and bottom

The north and south walls are **insulating** for heat and salt, and the bottom is
insulating and free-slip. (A prescribed depth-dipole heat flux on the N/S walls,
inherited from an earlier reference run, was carried for a long time; its net
contribution was $+0.0003$ PW and its per-cell tendency ~4000× smaller than the SST
restoring, so it was removed on 2026-10-21.)

---

## 5. Body forcings

### 5.1 Momentum sponges and deep Rayleigh drag

A linear (Rayleigh) damping is applied to the horizontal velocities,

$$
\mathbf{F}_{u,v} = -\big[\lambda_0\,m_{\rm sp}(x,y) + \lambda_d\,m_d(z)\big]\,(u,v).
$$

The boundary sponge $m_{\rm sp}$ damps velocities within 500 km of the east/west
and north/south walls (quadratic ramps), absorbing boundary reflections:

$$
m_{\rm sp}(x,y) = \alpha_{\rm EW}\max(w_W^2,w_E^2) + \alpha_{\rm NS}\max(w_S^2,w_N^2).
$$

The deep Rayleigh drag $m_d(z)$ is zero in the upper ocean (where the EUC lives)
and ramps to full strength below ≈ 550 m,

$$
m_d(z) = \mathrm{clamp}\!\Big(\frac{z_{\rm top}-z}{z_{\rm top}-z_{\rm bot}},\,0,\,1\Big),
\qquad z_{\rm top}=-300\,\mathrm{m},\ z_{\rm bot}=-550\,\mathrm{m}.
$$

While unphysical, this drag keeps the EUC from extending below ~300 m and removes
the zonal momentum that would otherwise accumulate at depth in a closed basin.

| $\lambda_0$ | sponge width | $\lambda_d$ | $z_{\rm top}$ | $z_{\rm bot}$ |
|---|---|---|---|---|
| $1/5$ days | 500 km (each wall) | $1/20$ days | −300 m | −550 m |

![Momentum sponge and deep drag](figures/05_momentum_sponge_deep_drag.png)

### 5.2 Poleward heat sink (north/south temperature relaxation)

In a closed basin the surface heat taken up in the cold tongue (~2.6 PW here) has
nowhere to go and the ocean warms indefinitely. To mimic the real
tropical→extratropical heat export, $T$ is relaxed toward a **cool extratropical
reference profile** within the north/south bands, **in the upper ocean only**:

$$
F_T^{\rm sink} = -\lambda_{\rm sink}\,m_{\rm NS}(y)\,m_z(z)\,\big(T - T_{\rm ref}(z)\big),
$$
$$
T_{\rm ref}(z) = T_{\rm deep}^{\rm ref}
+ (T_{\rm surf}^{\rm ref}-T_{\rm deep}^{\rm ref})\,
\tfrac12\Big[1 + \tanh\!\Big(\frac{z - z_{\rm ref}}{d_{\rm ref}}\Big)\Big],
\qquad
m_z(z) = \tfrac12\Big[1 + \tanh\!\Big(\frac{z - z_m}{d_m}\Big)\Big].
$$

Warm water carried to the subtropical edges is cooled there, exporting the tropical
heat uptake so the basin equilibrates. The depth mask $m_z$ (≈1 above 150 m, 0.5 at
300 m, ≈0 below 450 m) confines the sink to the upper ocean: without it the sink
removed a further 6.6 PW from 300–900 m band water — the water that feeds the
EUC — chilling the cold tongue to 14 °C and cooling the basin.

| $\lambda_{\rm sink}$ | $T_{\rm surf}^{\rm ref}$ | $T_{\rm deep}^{\rm ref}$ | $z_{\rm ref}$ | $d_{\rm ref}$ | $z_m$ | $d_m$ |
|---|---|---|---|---|---|---|
| $1/150$ days | 22 °C | 3.5 °C | −150 m | 120 m | −300 m | 75 m |

![Poleward heat sink](figures/06_poleward_heat_sink.png)

In the equilibrated spin-up the sink removes 2.95 PW (1.6 above 150 m, 1.4 below),
balancing the 2.6 PW surface restoring input to within the ~0.3 PW residual.

---

## 6. Subgrid closure and advection

Three closures are stacked; their diffusivities add.

**Convective adjustment with a small background.** A tiny background diffusivity
(to preserve the thermocline) plus a large diffusivity that fires only where the
column is statically unstable ($N^2<0$). A hydrostatic model cannot represent
convection, so this is the safety valve that removes static instabilities:

$$
\kappa_z = \begin{cases}\kappa_{\rm conv} & N^2 < 0\\[2pt]\kappa_{\rm bg} & N^2 \ge 0\end{cases},
\qquad \kappa_{\rm conv}=1.0,\ \kappa_{\rm bg}=10^{-5},\ \nu_{\rm bg}=10^{-3}\ \mathrm{m^2\,s^{-1}},\ \nu_{\rm conv}=0 .
$$

**Richardson-number-dependent shear mixing** (Oceananigans
`RiBasedVerticalDiffusivity`, the Pacanowski–Philander class). Both viscosity and
diffusivity are a ceiling value times a taper in the gradient Richardson number
$Ri = N^2/(\partial_z u^2 + \partial_z v^2)$:

$$
\nu = \nu_0\,\tau(Ri),\quad \kappa = \kappa_0\,\tau(Ri),\qquad
\tau(Ri) = \tfrac12\Big[1 - \tanh\!\Big(\frac{Ri - Ri_0}{Ri^{\delta}}\Big)\Big],
$$

with $\nu_0 = \kappa_0 = 10^{-2}\ \mathrm{m^2\,s^{-1}}$ (Prandtl number 1),
$Ri_0 = 0.1$, $Ri^{\delta} = 0.4$: the taper is 0.62 at $Ri=0$, 0.32 at $Ri=0.25$
and 0.01 at $Ri=1$. $Ri$ is smoothed with a five-point horizontal filter before the
taper is evaluated, which breaks the noisy-shear → noisy-mixing feedback that
produced grid-scale lines at coarse resolution. The scheme's convective term is
switched off ($\kappa^{\rm ca}=0$) because the convective adjustment above already
handles $N^2<0$. In the spun-up state the SEC–EUC shear layer sits at
$Ri\approx0.2$–0.3, where this gives $\nu\approx 4\times10^{-3}$. (A fit of this
form to the PP81 power law was considered and set aside: PP81 is 3–4× less
diffusive in that $Ri$ range; the plain Oceananigans values with $Pr=1$ were kept
for simplicity. Reducing $\kappa_0$ or $\nu_0$ alone did not change the SEC in
branch tests — the $Ri$ feedback compensates.)

**Horizontal mixing** is a **biharmonic viscosity on momentum only**,
$\partial_t\mathbf{u}_h \supset -\nu_4\,\nabla_h^4\mathbf{u}_h$ with
$\nu_4 = 10^{12}\ \mathrm{m^4\,s^{-1}}$: an e-folding of about an hour at
$2\Delta x$ and about 20 years at 1000 km. It is *not* applied to tracers, because
biharmonic diffusion is not sign-preserving and can drive salinity negative.

All vertical mixing is **implicit** in time, so the large convective diffusivity
does not constrain the time step.

**Advection** uses WENO schemes: order 3 for momentum and **order 5 for tracers**
(order 7 left the advected cold-tongue front ringing with $2\Delta x$ noise at
coarse resolution; order 5 adds enough dissipation to damp it while staying far less
diffusive than order 3).

---

## 7. Initial conditions

Both tracers start horizontally uniform with idealised tropical stratification: a
warm mixed layer over a sharp thermocline (temperature), and a subsurface salinity
maximum (halocline). Velocities start from rest.

$$
T_i(z) = T_d + (T_m-T_d)\,\tfrac12\!\left[1+\tanh\tfrac{z-z_{dt}}{\delta_{dt}}\right]
+ (T_s-T_m)\,\tfrac12\!\left[1+\tanh\tfrac{z-z_{th}}{\delta_{th}}\right]
$$

with $T_s=28$, $T_m=12$, $T_d=3.5$ °C, thermocline at $z_{th}=-140$ m
($\delta_{th}=55$ m), deep transition at $z_{dt}=-900$ m ($\delta_{dt}=450$ m).
The salinity profile $S_i(z)$ has an analogous two-tanh form ($S=34.6$ surface,
35.1 subsurface, 34.8 deep).

![Initial conditions](figures/07_initial_conditions.png)

---

## 8. Time stepping, output and checkpointing

The simulation uses an adaptive time step (CFL = 0.2, capped at 600 s; the cap is
active throughout). On one V100 the model runs ≈ 660 model days per wall-hour, so
a five-year spin-up takes ≈ 3 h and a one-year branch ≈ 35 min.

**Output** (all JLD2):

| File | Contents | Cadence |
|---|---|---|
| `toy_model.jld2` | 3-D $T,S,u,v,w,\partial_z T$, **time-averaged** (`AveragedTimeInterval`) | 30-day means (spin-up), 10-day means (branches) |
| `toy_model_surface.jld2` | top-cell $T,S,u,v,w,\partial_z T$ | instantaneous, 2 days |
| `toy_model_140W.jld2`, `_110W.jld2` | meridional $y$–$z$ sections of $T,S,u,v,w$ at $x = 0.625L_x$, $0.81L_x$ | instantaneous, 2 days |
| `toy_model_ymean.jld2` | $y$-mean $T$ and $u$ ($x$–$z$) | instantaneous, 2 days |
| `toy_model_momentum.jld2` | volume-integrated zonal momentum (upper/deep), $\int|u|$, $\int u^2$ below 300 m | daily |
| `checkpoint_iteration<N>.jld2` | full model state (`Checkpointer`), latest kept | every 100 days |

The time-averaged 3-D file is what the equilibrium diagnostics use; the 2-day
slices carry the high-frequency information (tropical instability waves, EUC
variability) at a small fraction of the size.

**Checkpoint / branch workflow.** `run_model` takes a `pickup` keyword: `false`
starts from the initial conditions; a checkpoint path restores the full state and
clock, after which `stop_time` is absolute (checkpoint time + branch length). A
pickup is a bit-identical continuation of the original run. Forcing parameters may
differ between the checkpointed run and the branch (that is how perturbation
experiments are made); grid and closures may not.

The run script `run_toy_model.jl` selects the configuration from the `BRANCH`
environment variable at submit time (`qsub -v BRANCH=…`), so several branches can
queue together; `BRANCH_YEARS` (default 1) sets every branch's length.

**Diagnostics.** Two families are computed or post-processed:

- **Momentum spin-down**: the volume-integrated zonal momentum split at 300 m, and
  $\int|u|\,\mathrm{d}V$ below 300 m. Both flattening confirms the deep ocean has
  spun down rather than accumulating momentum.
- **Heat budget**: the restoring flux is reconstructed from the saved SST
  (restoring in / out / net) and compared with the heat-content tendency
  $\mathrm{d}/\mathrm{d}t\int\rho c_p T\,\mathrm{d}V$; the gap is the sink. With
  the sink computed from the saved $T$ the budget closes to ±0.05 PW, which is also
  the check that the model conserves heat (§12).

---

## 9. Spin-up and equilibration

The model starts from rest with horizontally uniform stratification (§7). The
trade-wind stress is ramped on over the first 10 days and is steady thereafter.
Two timescales govern the approach to a statistically steady state.

**Dynamical adjustment — fast, with an overshoot.** Once the winds are up, the
surface flow and the equatorial wave field redistribute mass: warm water piles up
in the west, the thermocline tilts, and the subsurface zonal pressure gradient
drives the EUC. The adjustment is not monotonic: around days 100–200 the EUC
transiently shoals and balloons, nearly flipping the surface eastward, before the
wave field settles and the SEC re-establishes on top by day ~250. The momentum
integrals are at their equilibrium values by day ~400.

**Thermal adjustment — slow.** The heat content equilibrates over years, because
the SST restoring and the poleward sink must come into balance through the
thermocline. In the production spin-up the heat-content tendency falls to
$-0.3$ PW by day ~600 and stays there to the end of the five years (a residual of
$-7\ \mathrm{W\,m^{-2}}$ basin-average, accepted); the thermocline depth drifts by
~2 m per 100 days and SST by < 0.1 K per 100 days in the final year.

Production runs are therefore spun up for **1825 days (5 years)**, with the last
checkpoint at **day 1800**, and all experiments branch from that checkpoint.
Time-mean diagnostics are formed over the **last 100 days**; judging the EUC/SEC
from an early snapshot is misleading.

---

## 10. Resulting state

Production spin-up `runs/toy_model_spinup_2026-10-21_140641`, days 1710–1800 (last three
30-day means) for the figure; days 1725–1825 from the 2-day slices for the table.

| Quantity | Model | Observed (140°W) |
|---|---|---|
| EUC core at 140°W ($\|y\|\le100$ km) | 0.96 m s⁻¹ at 110 m | 1.0–1.2 m s⁻¹ at 110–130 m |
| EUC core at 110°W | 0.83 m s⁻¹ at 90 m | ~0.8 m s⁻¹ at ~90 m |
| Equatorial-band ($\pm250$ km) EUC | 0.53 ± 0.04 m s⁻¹ at 116 m | — |
| Surface zonal velocity on the equator | −0.13 m s⁻¹ (band ±250 km: −0.31) | −0.2 to −0.5 m s⁻¹ |
| Westward layer thickness | ~20 m | 40–60 m |
| 20 °C isotherm: west / 140°W / east | 211 / 153 / 50 m | ~170 / 120 / 50 m |
| SST: warm pool / east | 29.5 / 19.2 °C | 29 / 23 °C |
| Steric height W − E (equator) | 0.68 m | ~0.55 m |
| Basin-mean SSS | 35.0 | — |

![Equatorial cross section](figures/08_result_equatorial_section.png)

*Equatorial-band ($|y|\le250$ km) mean over the last 100 days of the spin-up.
Top: temperature with the 20 °C isotherm (white). Bottom: zonal velocity, the
eastward EUC beneath the westward SEC, and weak flow below 300 m.*

The EUC's core speed, depth, width (≈ 240 km half-width), eastward deceleration
and intraseasonal variability (±0.05 m s⁻¹ from tropical instability waves) are
realistic. Two known idealisations remain: the cold tongue is ~5 K colder than its
target (eastern-wall upwelling in a closed box outruns the restoring), and the
westward surface layer on the equator is thinner than observed (the EUC's upper
flank is lifted by the 2–4 m day⁻¹ equatorial upwelling). Eastward
counter-currents appear along both walls at ±12–14°, driven by the wind-stress curl
between the trade maximum and the weak-wind walls.

**Scaling check.** Pedlosky's (1987) inertial scales with the model's parameters
($\tau/\rho = 6.9\times10^{-5}\ \mathrm{m^2\,s^{-2}}$, reduced gravity
$\gamma_2 \approx 0.032\ \mathrm{m\,s^{-2}}$, $L = L_x$) give a thermocline depth
$H = (\tau L/\gamma_2)^{1/2} = 180$ m (model 150–210 m) and an inviscid jet speed
$U = (\gamma_2 L\tau)^{1/4} = 2.4$ m s⁻¹, of which the model's viscous EUC
realises ~40 %.

---

## 11. Experimental design: warming branches

All experiments are **branches from the spin-up checkpoint (day 1800)**, run
alongside a control branch of the same length; responses are perturbation minus
control over the same window, so any residual drift cancels. With the wind fixed,
"surface warming" is imposed through the SST restoring target (§4.1) and, in one
protocol, through the water subducted in the sink bands (§5.2). The run script
exposes three keywords: `ΔT_target` (added to both targets), `ΔT_east` (extra on
the cold-tongue target) and `ΔT_sink` (added to $T_{\rm surf}^{\rm ref}$).

| Branch | Protocol | $T_{\rm warm}$ | $T_{\rm cold}$ | $T_{\rm surf}^{\rm ref}$ | Isolates |
|---|---|---|---|---|---|
| `control` | — | 29.5 | 24.0 | 22 | reference |
| `warmA1` / `warmA2` / `warmA4` | **A** surface warming | +1 / +2 / +4 | +1 / +2 / +4 | unchanged | warmer surface over unchanged EUC source water → stratification ↑ |
| `warmB2` | **B** ventilated warming | +2 | +2 | +2 | the subducted source water warms too → smaller stratification change; A − B is the role of source-water warming |
| `warmC2` | **C** pattern warming | +2 | +4 | unchanged | zonal SST gradient weakened (CMIP-like eastern-enhanced warming) → zonal pressure gradient |

![Warming experiment forcing profiles](figures/09_warming_experiments.png)

*Left: the SST restoring target along the equator for each branch. Middle: the
imposed change relative to the control (B2 coincides with A2 at the surface).
Right: the sink's reference profile, unchanged in A and C and warmed by 2 K at the
surface in B2.*

**Why the target and not a flux.** With $\lambda_T \approx 164\ \mathrm{W\,m^{-2}\,K^{-1}}$
an added flux of 20 W m⁻² would move SST by 0.1 K; in a restoring framework the
warming must be imposed through the target.

**Analysis plan.** For each branch: thermocline depth ($z_{20}$ or the $N^2$
maximum), EUC core depth, speed, half-width and transport, the temperature
(isopycnal) of the core, and the reduced gravity $\gamma_2$ diagnosed from the
run's own thermocline. The inertial expectation (Pedlosky 1987, §10) is
$H \propto \gamma_2^{-1/2}$ and $U \propto \gamma_2^{1/4}$: for +2 K on a 13 K
thermocline contrast, a ~8 % shoaling and ~4 % strengthening. The A1/A2/A4 series
tests the exponents on log–log axes; tracking the core isotherm separates an
adiabatic response (the core rides the thermocline up) from a diabatic one (the
core's temperature changes). One-year branches give the transient response; five-year
branches (`BRANCH_YEARS=5`, with a five-year control) the adjusted one.

---

## 12. Implementation notes

- **Oceananigans versions.** Gadi runs Oceananigans 0.96.35 (CUDA.jl 5.8, V100);
  the laptop has 0.107. In 0.96.35 a *continuous-form* tracer forcing with
  `field_dependencies = :T` receives the wrong field (≈ 0), which silently turned
  the poleward sink into a constant +15.6 PW heater in the first GPU runs. All body
  forcings (sink, sponges) are therefore written in **discrete form**, indexing
  `model_fields` by name, which is immune to the problem and bit-identical to the
  continuous form where the latter works. Boundary conditions with
  `field_dependencies` (SST and salinity restoring) are unaffected.
- **Heat-conservation check.** The closed budget (§8) is the regression test:
  $\mathrm{d}HC/\mathrm{d}t$ minus (surface restoring + sink) should be within
  ±0.1 PW once the first month's adjustment has passed.
- **Reading output off Gadi.** The JLD2 files from 0.96.35 do not deserialise with
  `FieldTimeSeries` under 0.107; read `timeseries/<var>/<iteration>` as raw arrays
  (4-cell halos, `Float32`) and rebuild the coordinates from §1.

---

### Parameter summary

| Group | Parameters |
|---|---|
| Grid | $L_x=15{,}000$ km, $L_y=\pm15°$, $L_z=2000$ m; $600\times150\times60$ (25 × 22 km; 9.5 m upper 300 m); GPU |
| Rotation | β-plane, $\beta=2.29\times10^{-11}$ m⁻¹s⁻¹ |
| Wind | $u_{10}$: 7 m s⁻¹ equator, 8 m s⁻¹ bump at ±10° (peak ≈ 7.7 at ±7°), quartic wall drop 4 m s⁻¹ (wall ≈ 3.6); $C_D=1.2\times10^{-3}$, $\rho_a=1.2$; $\tau_x$ = 0.071 N m⁻² at the equator; 10-day ramp; resting-ocean stress |
| Heat flux | pure SST restoring, $\lambda_T=4\times10^{-5}$ m s⁻¹; $T_{\rm target}$: 29.5 → 24 °C, front at 0.6 $L_x$, width 0.18 $L_x$ |
| Salinity | uniform evaporation $1\times10^{-8}$ m s⁻¹ + restoring $5\times10^{-7}$ m s⁻¹ to 34.6 g kg⁻¹ |
| Walls, bottom | insulating; free-slip bottom |
| Momentum drag | sponge $1/5$ day within 500 km of all walls; deep drag $1/20$ day ramping over 300–550 m |
| Heat sink | $1/150$ day toward $T_{\rm ref}(z)$ (22 °C surface, 3.5 °C deep, $z_{\rm ref}=-150$ m, $d_{\rm ref}=120$ m) in the N/S 500-km bands; depth mask 0.5 at 300 m, width 75 m |
| Closure | convective adjustment ($\kappa_{\rm conv}=1$, $\kappa_{\rm bg}=10^{-5}$, $\nu_{\rm bg}=10^{-3}$) + Ri scheme ($\nu_0=\kappa_0=10^{-2}$, $Ri_0=0.1$, $Ri^\delta=0.4$, 5-point filter) + biharmonic $\nu_4=10^{12}$ m⁴s⁻¹ (momentum) |
| Advection | WENO-3 momentum, WENO-5 tracers |
| Time step | CFL 0.2, max 600 s |
| Output | 3-D means 30 d (spin-up) / 10 d (branches); 2-day slices; checkpoints every 100 d |
| Spin-up | 1825 d from rest; branches from the day-1800 checkpoint |
| Experiments | A: targets +1/+2/+4 K; B: A2 + sink ref +2 K; C: warm pool +2, cold tongue +4 K |
