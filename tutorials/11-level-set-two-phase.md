# Tutorial 11 — Two-Phase Flow with Level Sets

**Goal:** Make an interface move. Track it, keep it clean, apply surface tension
as a force, and understand the two failure modes (mass loss and spurious
currents) that every two-phase code fights.

**Code:** `code/lib/level_set.m`, `code/lib/curvature.m`, `code/ch11/static_droplet.m`, `code/ch11/rising_bubble.m`

---

## 11.1 How to represent a moving interface

| Method | Idea | Mass conservation | Curvature | Topology changes |
|---|---|---|---|---|
| **Front tracking** | explicit marker points | excellent | excellent | manual, painful |
| **VOF** | volume fraction per cell | **exact** | poor (needs reconstruction) | automatic |
| **Level set** | signed distance $\phi$ | poor (leaks) | **excellent** | automatic |
| **Phase field** | diffuse order parameter | good | good | automatic |

This course uses **level sets**, because curvature accuracy is what surface
tension needs (Tutorial 10.6) and because breakup/coalescence must be automatic
for Tutorial 12's droplet generation.

The price is mass conservation, and it is a real price. We measure it rather
than hope.

*(Production codes often use **CLSVOF** — VOF for the volume, level set for the
curvature — getting both. Worth knowing that's the state of the art.)*

---

## 11.2 The level-set idea

Define $\phi(\mathbf{x},t)$ as the **signed distance** to the interface:

$$\phi > 0 \text{ in fluid 1}, \qquad \phi < 0 \text{ in fluid 2}, \qquad \phi = 0 \text{ on the interface}$$

The interface is the zero contour — never stored explicitly, so it can merge and
split without any bookkeeping.

Signed distance means $|\nabla\phi| = 1$, which is what makes everything
downstream clean:

$$\mathbf{n} = \frac{\nabla\phi}{|\nabla\phi|} = \nabla\phi, \qquad \kappa = \nabla\cdot\mathbf{n} = \nabla^2\phi$$

Normals and curvature become trivial — *provided* $\phi$ stays a distance
function, which is the whole difficulty (§11.3).

Advect it with the flow:

$$\frac{\partial\phi}{\partial t} + \mathbf{u}\cdot\nabla\phi = 0$$

Use a **high-order scheme**. $\phi$ is smooth, so the wiggle problem of Tutorial 5
is absent, but numerical diffusion here directly *moves the interface*. We use
WENO or the van Leer scheme from Tutorial 5.

### Fluid properties from $\phi$

Density and viscosity jump across the interface. A sharp jump makes the pressure
solve badly conditioned, so smear it over ~1.5 cells with a smoothed Heaviside:

$$H_\varepsilon(\phi) = \begin{cases} 0 & \phi < -\varepsilon\\ \frac{1}{2}\left(1 + \frac{\phi}{\varepsilon} + \frac{1}{\pi}\sin\frac{\pi\phi}{\varepsilon}\right) & |\phi|\le\varepsilon\\ 1 & \phi>\varepsilon\end{cases}$$

$$\rho = \rho_2 + (\rho_1-\rho_2)H_\varepsilon(\phi), \qquad \mu = \mu_2 + (\mu_1-\mu_2)H_\varepsilon(\phi)$$

with $\varepsilon = 1.5\Delta x$. Smaller $\varepsilon$ is sharper but noisier;
larger is smoother but smears the physics. 1.5 is the standard compromise.

---

## 11.3 Reinitialization: the part that makes or breaks it

**The problem.** Advection does not preserve $|\nabla\phi|=1$. In a shearing
flow, $\phi$ stretches in some places and compresses in others. After a few
hundred steps it is no longer a distance function, and then:

- $\kappa = \nabla^2\phi$ is wrong (it's only valid when $|\nabla\phi|=1$)
- the smeared interface width varies unpredictably
- surface tension forces become garbage

**The fix.** Periodically restore the distance property by marching this to
steady state in pseudo-time $\tau$:

$$\frac{\partial\phi}{\partial\tau} = \text{sign}(\phi_0)\left(1 - |\nabla\phi|\right)$$

This drives $|\nabla\phi|\to1$ while holding the zero contour fixed — because
$\text{sign}(0)=0$, the interface itself doesn't move. Or shouldn't; in practice
the *discrete* sign function isn't exactly zero at the interface and the contour
drifts slightly. That drift is the main source of level-set mass loss.

Use the smoothed sign function

$$\text{sign}_\varepsilon(\phi_0) = \frac{\phi_0}{\sqrt{\phi_0^2 + |\nabla\phi_0|^2\Delta x^2}}$$

which reduces (but does not eliminate) the drift.

**Practical settings:**
- Reinitialize every 20 advection steps (this course's default).
- 3–5 pseudo-time iterations each time, with $\Delta\tau = 0.5\Delta x$.
- Only near the interface (a band of ±5 cells) — far away nobody cares, and it's
  cheaper.

**Too much reinitialization loses mass; too little corrupts curvature.** There is
no universally right setting. Tune it by watching the area diagnostic.

### How bad is the drift? Bad enough to destroy the droplet

Measured on a circle: each reinitialization call moves the enclosed area by
about **0.1%**. That sounds harmless. It is not. A static-droplet run of 30,000
steps reinitializing every 10 steps makes 3,000 calls, and $0.999^{3000}\approx0$
— the droplet **completely dissolves**, with the area error reaching $-1.0$
(total loss). This actually happened while building this course, and it looks
exactly like a physics result until you check the diagnostic.

So this course applies a **global area correction** after each reinitialization:
shift $\phi$ by a constant $\delta$ chosen so the enclosed area returns to
target. Since $dA/d\delta = \oint dl$ (the perimeter),

$$\delta = \frac{A_{target} - A}{\text{perimeter}}$$

applied as two or three Newton steps (`level_set('fixarea',...)`). With it, the
measured area error over the same run drops from $-1.0$ to $-1.5\times10^{-6}$.

**Two honest caveats.** It fixes *total* volume, not its local distribution, and
it moves the whole interface uniformly rather than where the loss occurred — a
patch on a known weakness, not a cure. And it must be switched **off** when the
domain has inflow/outflow and the dispersed volume is genuinely supposed to grow
(droplet generation, Tutorial 12) — `opt.conserveMass = false` there. Leaving it
on would fight the injection and quietly produce a completely wrong droplet size.

---

## 11.4 Surface tension: the CSF model

Surface tension acts *on* the interface, but our grid has no interface — just a
smeared region. The **Continuum Surface Force** model (Brackbill 1992) converts
the surface force into an equivalent volumetric force:

$$\mathbf{f}_{st} = \sigma\kappa\,\delta_\varepsilon(\phi)\,\mathbf{n} = \sigma\kappa\,\delta_\varepsilon(\phi)\nabla\phi$$

with the smoothed delta function $\delta_\varepsilon = dH_\varepsilon/d\phi$:

$$\delta_\varepsilon(\phi) = \frac{1}{2\varepsilon}\left(1 + \cos\frac{\pi\phi}{\varepsilon}\right) \text{ for } |\phi|\le\varepsilon,\quad 0 \text{ otherwise}$$

Add $\mathbf{f}_{st}/\rho$ to the momentum predictor. That's the whole coupling.

### Spurious currents — the characteristic failure

Take a **static droplet** in zero gravity with no flow. Exactly, nothing should
happen: surface tension is balanced by the pressure jump, and $\mathbf{u}=0$
forever.

Numerically, you get small vortices at the interface that never decay. These
**spurious** (or parasitic) currents arise because the discrete
$\sigma\kappa\delta\nabla\phi$ and the discrete $\nabla p$ do not cancel exactly
— they're computed by different stencils at different locations.

Magnitude scales as

$$|u_{spurious}| \sim \frac{\sigma}{\mu}\times(\text{a number depending on your discretization})$$

**This is the standard benchmark for a two-phase code**, and `static_droplet.m`
runs it. What to check:

- Spurious velocity should be **small compared to your physical velocity**. Report
  the ratio, don't just eyeball the plot.
- It should **not grow** with time.
- Refining the grid should reduce it (slowly — this is a hard problem).
- The pressure jump across the interface should match $\sigma/R$ (2D!).

If your spurious currents are comparable to the flow you're trying to simulate,
your droplet results are noise. The cure is a **balanced-force** discretization:
compute $\kappa$ and $\nabla\phi$ at the *same* locations where $\nabla p$ lives
(i.e. on the faces), so the two discrete gradients can cancel. `ns2phase.m`
does this, and `static_droplet.m` quantifies how well.

### Measured results — including the uncomfortable one

Actual output from `static_droplet.m` (water-in-oil, $\sigma=5$ mN/m,
$\mu_c=10$ mPa·s, $R=40$ µm):

```
   N   cells/R   |u|/(sigma/mu)   dp err     area err
  32      6.4         1.83e-02      2.6%    -3.08e-06
  64     12.8         1.16e-02      1.0%    -1.50e-06
```

The pressure jump converges toward $\sigma/R$ (2.6% → 1.0%) and mass is
conserved to $10^{-6}$. Those pass.

**The spurious currents do not.** At $N=64$, $1.16\times10^{-2}\times(\sigma/\mu)
= 5.8\times10^{-3}$ m/s — which is **six times larger** than the 1 mm/s flow a
microchannel simulation would be trying to resolve. On this grid, with these
fluid properties, computed droplet trajectories and breakup times would be
dominated by discretization error rather than physics.

That is a real limitation of this implementation, stated plainly because it
determines how you should read Tutorial 12: **trust the regimes and the trends,
not the absolute droplet volumes.** Your options are to refine (converging, but
slowly), to choose fluids with a slower capillary scale $\sigma/\mu$ (a more
viscous oil or more surfactant — often physically legitimate), or to move to a
sharper curvature method such as height functions or CLSVOF, which is the real
fix and is beyond this course.

`static_droplet.m` prints this verdict explicitly rather than leaving you to do
the arithmetic. A spurious-current number is meaningless in isolation; it only
means something relative to the flow you intend to simulate.

---

## 11.5 The capillary timestep — your new smallest number

Surface tension supports capillary waves. Resolving them explicitly requires
(Brackbill et al. 1992):

$$\boxed{\Delta t \le \sqrt{\frac{(\rho_1+\rho_2)\Delta x^3}{4\pi\sigma}}}$$

Note the scaling: $\Delta t_{cap}\propto\Delta x^{3/2}$ while
$\Delta t_{visc}\propto\Delta x^2$. So as you refine, the *viscous* limit
eventually wins — but the capillary limit has a brutal constant at the
microscale, and at practical grid sizes it is the one that binds.

Water/air, $\Delta x = 1$ µm:

$$\Delta t \le \sqrt{\frac{1000\times10^{-18}}{4\pi\times0.072}} = 3.3\times10^{-8}\ \text{s}$$

Compare Tutorial 7's viscous limit of $2.5\times10^{-7}$ s. **Surface tension is
~8× more restrictive at this grid size**, and it is now your binding constraint.
To simulate 0.1 s of droplet formation you need ~3 million steps.

(Because of the differing exponents, the crossover is at
$\Delta x \approx \left[\frac{4\pi\sigma}{(\rho_1+\rho_2)}\cdot\frac{1}{(2\nu)^2}(\tfrac{1}{2})^{2}\right]^{-1}$
— rather than memorize that, just compute both limits in code and let `min` pick.
`ns2phase.m` prints which one is active, which is more useful than knowing the
formula.)

This is why droplet simulations are expensive, and why people use implicit or
semi-implicit surface tension for production work. For learning, accept the cost
and run small domains.

Your timestep is now the minimum of **four** limits:

```matlab
dt = safety * min([ dx/max(abs(u(:))), ...                    % CFL
                    0.5/nu/(1/dx^2+1/dy^2), ...               % viscous
                    sqrt((rho1+rho2)*dx^3/(4*pi*sigma)), ...  % capillary
                    dtMax ]);
```

---

## 11.6 Validation ladder

Do these **in order**. Each one isolates a different piece, and skipping ahead
means debugging several things at once.

**1. Advection only** (`level_set.m` self-test). Zalesak's disk: a slotted disk
in solid-body rotation. After one full revolution it should return to its
initial shape. Measures interface-tracking error alone — no surface tension, no
coupling.

**2. Curvature** (`young_laplace.m`, Tutorial 10). A circle has $\kappa=1/R$
exactly. Check the convergence rate.

**3. Static droplet** (`static_droplet.m`). No flow, no gravity. Checks:
- pressure jump $= \sigma/R$ (2D)
- spurious currents small and not growing
- mass conserved

**4. Rising bubble** (`rising_bubble.m`). Buoyancy vs drag vs surface tension.
Compare terminal velocity against Stokes' law for a bubble:
$$U_t = \frac{2}{3}\frac{\Delta\rho\, g R^2}{\mu}\ \ \text{(clean bubble)}, \qquad U_t = \frac{2}{9}\frac{\Delta\rho\, g R^2}{\mu}\ \ \text{(rigid sphere)}$$
The factor differs because a clean bubble has a mobile interface. Which one you
get is itself a diagnostic of whether your interface is behaving.

**5. Droplet in shear.** Taylor deformation parameter
$D = (L-B)/(L+B) \approx \frac{19\lambda+16}{16\lambda+16}Ca$ for small $Ca$,
where $\lambda = \mu_{in}/\mu_{out}$.

Only after all five should you trust Tutorial 12's droplet generation.

---

## 11.7 Honest limitations of this implementation

State these when you report results:

- **Mass loss.** Measured and reported every run. Typically 1–5% over a long
  simulation on a moderate grid. Unacceptable for precise droplet volumes;
  use VOF or CLSVOF if that's your quantity of interest.
- **Spurious currents.** Reduced by balanced-force discretization but not
  eliminated. Report the ratio to physical velocity.
- **Contact lines.** Not implemented here. A moving contact line on a no-slip
  wall is a genuine *singularity* in continuum mechanics (infinite viscous
  dissipation), and resolving it requires a slip model with a microscopic
  length. This is a research topic, not an oversight to patch. For droplet
  generation in a wetted channel it doesn't matter — the continuous phase coats
  the walls and the interface never touches them.
- **2D.** A 2D "droplet" is an infinite cylinder. Curvature is $1/R$, not $2/R$;
  breakup dynamics differ; and the Rayleigh–Plateau instability that drives real
  jet breakup is fundamentally 3D. 2D gets the qualitative regimes right and the
  quantitative numbers wrong.
- **Density ratio.** Water/air is 1000:1, which is numerically hard. Water/oil
  (~1:1) is much easier and is what most droplet microfluidics uses anyway.
  Start there.

---

## Exercises

**11.1** Run the Zalesak disk test in `level_set.m`. Measure the shape error
after one revolution. Then disable reinitialization and repeat — how much worse?

**11.2** Run `static_droplet`. Verify $\Delta p = \sigma/R$ (not $2\sigma/R$ —
check you understand why). Report the max spurious velocity as a fraction of
$\sigma/\mu$.

**11.3** Run `static_droplet` on 32², 64², 128². Plot spurious current magnitude
vs grid. What's the observed convergence rate? *(Expect it to be poor. This is
the honest state of the art, not a defect in this code.)*

**11.4** Track mass over 5000 steps of `static_droplet`. Plot the loss. Now vary
the reinitialization frequency (every 2, 5, 10, 20 steps) and find the best
setting for your case.

**11.5** Run `rising_bubble`. Compare terminal velocity to both Stokes formulas.
Which does it match? What does that tell you about the interface condition your
code is producing?

**11.6** Derive the capillary timestep limit yourself: the capillary wave
dispersion relation is $\omega^2 = \sigma k^3/(\rho_1+\rho_2)$; require that the
shortest resolvable wave ($k = \pi/\Delta x$) be resolved in time.

**11.7 (the one that matters)** Set up a droplet in simple shear at
$Ca = 0.05$, $\lambda = 1$. Measure the deformation parameter and compare to
Taylor's formula. This validates the *coupling* of flow, interface, and surface
tension all at once — and it's the last check before Tutorial 12.

---

**Next:** [Tutorial 12 — Droplet generation](12-droplet-generation.md).
Everything comes together: geometry, two-phase flow, surface tension, and a
design question with a real answer.
