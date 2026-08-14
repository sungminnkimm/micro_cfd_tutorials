# Appendix C — Notation, Symbols, and Formula Reference

---

## C.1 Symbols

| Symbol | Meaning | Units |
|---|---|---|
| $\mathbf{u} = (u,v)$ | velocity | m/s |
| $p$ | pressure | Pa |
| $\rho$ | density | kg/m³ |
| $\mu$ | dynamic viscosity | Pa·s |
| $\nu = \mu/\rho$ | kinematic viscosity | m²/s |
| $\sigma$ | surface / interfacial tension | N/m |
| $D$ | mass diffusivity | m²/s |
| $c$ | scalar concentration | — |
| $\phi$ | level-set function (T11) / pressure correction (T7) | m / m²/s |
| $\kappa$ | total curvature | 1/m |
| $\theta$ | contact angle | deg |
| $\psi$ | streamfunction | m²/s |
| $L, w, h$ | length, width, height | m |
| $\Delta x, \Delta y, \Delta t$ | grid spacing, timestep | m, s |

**Watch out:** $\phi$ means two different things in this course — the pressure
correction in Tutorial 7 and the level set in Tutorial 11. They never appear in
the same equation, but the collision is standard in the literature so it's worth
being alert to.

---

## C.2 Dimensionless groups

| Group | Definition | Compares | Micro value |
|---|---|---|---|
| $Re$ | $\rho U L/\mu$ | inertia / viscous | $10^{-3}$–10 |
| $Pe$ | $UL/D$ | advection / mass diffusion | 10–$10^4$ |
| $Sc$ | $\nu/D$ | momentum / mass diffusivity | ~2000 (water) |
| $Ca$ | $\mu U/\sigma$ | viscous / surface tension | $10^{-5}$–$10^{-2}$ |
| $Bo$ | $\rho g L^2/\sigma$ | gravity / surface tension | $10^{-3}$ |
| $We$ | $\rho U^2L/\sigma = Re\,Ca$ | inertia / surface tension | $10^{-6}$ |
| $Oh$ | $\mu/\sqrt{\rho\sigma L}$ | viscous / inertia-capillary | — |
| $Kn$ | $\lambda/L$ | mean free path / size | <$10^{-3}$ (liquids) |
| $De$ | $Re\sqrt{h/2R}$ | Dean (secondary flow in bends) | <1 |
| $Pe_{cell}$ | $\|u\|\Delta x/D$ | **numerical**: wiggle threshold is 2 | — |
| $C$ | $\|u\|\Delta t/\Delta x$ | **numerical**: Courant, CFL limit 1 | — |
| $r$ | $D\Delta t/\Delta x^2$ | **numerical**: diffusion number, limit 1/2 | — |

The last three are properties of your *discretization*, not your flow. Confusing
$Pe$ with $Pe_{cell}$ is a common and consequential slip.

---

## C.3 Analytic results used for validation

**Poiseuille, parallel plates, gap $h$, $G = -dp/dx$:**
$$u(y) = \frac{G}{2\mu}\left(\frac{h^2}{4}-y^2\right),\quad
\bar u = \frac{Gh^2}{12\mu},\quad \frac{u_{max}}{\bar u} = \frac{3}{2},\quad
Q = \frac{Gh^3}{12\mu}$$

**Rectangular duct correction** (multiply the 2D result): see
`duct_correction.m`. $w/h=1\to0.42$; $w/h=10\to0.94$.

**Entrance length:** $L_e \approx 0.06\,Re\,h$.

**Taylor–Aris**, plates separated by $2h$:
$$D_{eff} = D\left(1+\frac{2}{105}Pe_h^2\right),\quad Pe_h = \frac{\bar u h}{D}$$
(circular tube: $1/48$ instead of $2/105$). Valid for $t \gg h^2/(\pi^2 D)$.

**Diffusive mixing length:** $L_{mix}\approx 0.3\,Pe\,w$ for a half-seeded
channel to $M<0.05$ — the calibrated version of the $Pe\,w$ rule of thumb.

**Young–Laplace:** $\Delta p = \sigma\kappa$. Sphere $2\sigma/R$; **2D cylinder
$\sigma/R$** — this course is 2D.

**Capillary rise:** $h = 2\sigma\cos\theta/(\rho g r)$.
**Capillary length:** $\ell_c = \sqrt{\sigma/\rho g} \approx 2.7$ mm (water).
**Entry pressure, 2D channel width $w$:** $\Delta p = 2\sigma\cos\theta/w$.

**Stokes settling** (3D sphere): rigid $U=\frac{2}{9}\frac{\Delta\rho gR^2}{\mu}$;
clean bubble $\frac{2}{3}$; Hadamard–Rybczynski for ratio $\lambda$:
$U = \frac{2}{3}\frac{\Delta\rho gR^2}{\mu}\frac{\lambda+1}{3\lambda+2}$.

**Garstecki (T-junction, squeezing):** $L/w = 1+\alpha Q_d/Q_c$, $\alpha\approx1$.

**Taylor deformation (droplet in shear, small $Ca$):**
$D = \frac{19\lambda+16}{16\lambda+16}Ca$.

---

## C.4 Timestep limits (all four)

```matlab
dt_cfl   = dx / max(abs(u(:)));
dt_visc  = 0.5/nu / (1/dx^2 + 1/dy^2);            % 2-D
dt_diff  = 0.5/D  / (1/dx^2 + 1/dy^2);            % scalar
dt_cap   = sqrt((rho1+rho2)*dx^3/(4*pi*sigma));
dt = safety * min([dt_cfl dt_visc dt_diff dt_cap]);   % safety ~ 0.4-0.5
```

Which binds, at the microscale: **viscous** for single-phase momentum;
**capillary** once surface tension is present; **CFL** for high-$Pe$ scalars.

---

## C.5 Grid conventions used throughout

**Cell-centred, index order $(i,j) = (x,y)$, column-major flattening**
$k = i + (j-1)n_x$.

```matlab
[X,Y] = ndgrid(xc, yc);        % NOT meshgrid
contourf(xc, yc, F.');         % transpose to plot
```

**Staggered (MAC) sizes:**

| | size | on the walls? |
|---|---|---|
| `p` | `nx × ny` | no (ghosts everywhere) |
| `u` | `(nx+1) × ny` | **yes** at $x=0,L_x$; ghosts at $y$ walls |
| `v` | `nx × (ny+1)` | **yes** at $y=0,L_y$; ghosts at $x$ walls |

**Ghost-cell formulas:**

```
Dirichlet (value w):  ghost = 2*w - inner
Neumann (flux g):     ghost = inner -+ g*h     (sign from the outward normal)
No-slip wall:         ghost = -inner
Free slip:            ghost = +inner
```

---

## C.6 Fluid properties (20 °C)

| Fluid | $\rho$ [kg/m³] | $\mu$ [Pa·s] | $\sigma$ vs air [N/m] |
|---|---|---|---|
| Water | 998 | 1.00e-3 | 0.072 |
| Ethanol | 789 | 1.20e-3 | 0.022 |
| Glycerol | 1261 | 1.41 | 0.064 |
| Mineral oil | 850 | 3.0e-2 | 0.030 |
| Fluorinated oil (HFE-7500) | 1614 | 1.24e-3 | 0.016 |
| Air | 1.2 | 1.8e-5 | — |

**Interfacial tensions:** water/hexadecane 0.052; water/mineral oil + Span80
~0.005; water/fluorinated oil + surfactant ~0.005.

**Diffusivities in water:** small ion 2e-9; small molecule (fluorescein) 4.3e-10;
protein (BSA) 6e-11; 1 µm bead 4e-13 m²/s.

---

## C.7 Code map

| File | Purpose | Introduced |
|---|---|---|
| `lib/mac_grid.m` | staggered grid descriptor | T6 |
| `lib/ns_solver.m` | incompressible NS, projection method | T7 |
| `lib/ns_defaults.m` | solver options | T7 |
| `lib/ns2phase.m` | two-phase NS with CSF | T11 |
| `lib/ns2_defaults.m` | two-phase options | T11 |
| `lib/level_set.m` | interface tracking, reinit, area fix | T11 |
| `lib/curvature.m` | $\kappa$ and normals from $\phi$ | T10 |
| `lib/scalar_step.m` | 2D TVD advection–diffusion | T9 |
| `lib/advect_tvd.m` | 1D TVD advection | T5 |
| `lib/limiter.m` | flux limiters | T5 |
| `lib/make_mask.m` / `set_solid.m` | geometry masks | T8 |
| `lib/mixing_index.m` | mixing metric | T9 |
| `lib/duct_correction.m` | 2D → 3D flow-rate correction | T6 |
| `lib/streamfun.m` | streamfunction | T7 |
| `lib/padarray_neumann.m` | zero-flux ghost ring | T4 |
| `ch03/poisson2d.m`, `poisson_solve.m` | Poisson operator + solve | T3 |
| `ch02/diff_matrix.m` | 1D derivative matrices | T2 |
| `ch01/scaling.m` | dimensionless groups | T1 |
| `ch13/mms_driver.m`, `gci.m` | verification tools | T13 |
