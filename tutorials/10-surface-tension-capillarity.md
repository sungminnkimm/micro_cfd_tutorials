# Tutorial 10 — Surface Tension and Capillarity

**Goal:** Understand the force that dominates every interface at the microscale,
and learn to compute curvature — the quantity everything in Tutorials 11 and 12
depends on.

**Code:** `code/ch10/young_laplace.m`, `code/ch10/capillary_rise.m`, `code/ch10/meniscus.m`

---

## 10.1 Where surface tension comes from

A molecule in the bulk is pulled equally in all directions. A molecule at the
surface has neighbours on one side only, so it is pulled inward. Creating new
surface therefore costs energy, and $\sigma$ [J/m² = N/m] is that cost per unit
area. A surface behaves like a stretched elastic sheet.

Two readings of the same number, both useful:

- **Energy per area** [J/m²] — for equilibrium shapes and contact angles.
- **Force per length** [N/m] — for force balances on a contact line.

| Interface | $\sigma$ [mN/m] |
|---|---|
| Water / air | 72 |
| Water / air + surfactant | 25–40 |
| Water / hexadecane | 52 |
| Water / mineral oil + span80 | 5 |
| Ethanol / air | 22 |
| Mercury / air | 485 |

That surfactant row matters practically: adding surfactant drops $\sigma$ by
2–10×, which changes $Ca$ by the same factor and can move you between droplet
regimes (Tutorial 12). Surfactants are the main knob available for tuning
droplet generation, and they also stabilize droplets against coalescence.

---

## 10.2 Young–Laplace: curvature makes pressure

A curved interface supports a pressure jump:

$$\boxed{\Delta p = p_{in} - p_{out} = \sigma\left(\frac{1}{R_1} + \frac{1}{R_2}\right) = \sigma\kappa}$$

where $\kappa$ is the **total curvature** (sum of the two principal curvatures —
note this convention: some texts use *mean* curvature $H = \kappa/2$, giving
$\Delta p = 2\sigma H$; both are correct, and mixing them is a factor-of-2 bug).

Special cases you should know cold:

| Shape | $\kappa$ | $\Delta p$ |
|---|---|---|
| Sphere, radius $R$ | $2/R$ | $2\sigma/R$ |
| Cylinder, radius $R$ | $1/R$ | $\sigma/R$ |
| **2D circle** (our simulations) | $1/R$ | $\sigma/R$ |
| Flat | 0 | 0 |
| Saddle, $R_1 = -R_2$ | 0 | 0 |

**Note the third row.** Our 2D simulations model an infinite cylinder, so a
"droplet" of radius $R$ has $\Delta p = \sigma/R$, not $2\sigma/R$. Tutorial 11's
validation checks against $\sigma/R$ for exactly this reason, and comparing a 2D
result against the spherical formula is a very common factor-of-2 error.

### Run the numbers

A 50 µm water droplet in air:

$$\Delta p = \frac{2\times0.072}{25\times10^{-6}} = 5760\ \text{Pa} = 58\ \text{mbar}$$

That is *enormous* compared to the viscous pressures driving your flow — Tutorial
1 computed 0.67 Pa across an 800 µm channel. **Capillary pressure exceeds viscous
pressure by four orders of magnitude.** Practical consequences:

- A bubble stuck in a channel is nearly impossible to push out.
- Interfaces snap to their equilibrium shape essentially instantly.
- Droplets are spheres unless something strong deforms them.
- Numerically: surface tension is **stiff**, and it will set your timestep
  (Tutorial 11).

---

## 10.3 Contact angle and wetting

Where the interface meets a solid, the three surface energies balance
(**Young's equation**):

$$\sigma_{SV} = \sigma_{SL} + \sigma_{LV}\cos\theta \quad\Longrightarrow\quad \cos\theta = \frac{\sigma_{SV}-\sigma_{SL}}{\sigma_{LV}}$$

- $\theta < 90°$: **wetting** (hydrophilic). Water in glass or oxidized PDMS.
- $\theta > 90°$: **non-wetting** (hydrophobic). Water in native PDMS ($\theta\approx105°$).
- $\theta \to 0$: complete wetting, liquid spreads into a film.

**Surface chemistry is a design variable, not a material constant.** In droplet
microfluidics you *must* have the continuous phase wet the walls — otherwise the
dispersed phase touches the wall, sticks, and your droplet train fails. Making
water-in-oil droplets requires hydrophobic channels; oil-in-water requires
hydrophilic. PDMS is natively hydrophobic; oxygen-plasma treatment makes it
hydrophilic (temporarily — it recovers over hours to days, which is a real and
frequently-cursed experimental problem).

### Hysteresis

Real surfaces have a range of contact angles, not one. The interface advances at
$\theta_A$ and recedes at $\theta_R$, with $\theta_A > \theta_R$; the difference
(hysteresis) comes from roughness and chemical heterogeneity and can be 10–50°.

Most CFD, this course included, uses a single static angle. **This is a
genuine model limitation, not a numerical detail** — hysteresis is what pins
droplets on surfaces and what makes real contact-line motion history-dependent.
If contact-line behaviour is central to your problem, know that a static-angle
simulation will not capture it.

---

## 10.4 Capillary rise: the classic validation

Liquid in a tube of radius $r$ rises until capillary force balances weight:

$$\boxed{h = \frac{2\sigma\cos\theta}{\rho g r}}$$

Water in a 100 µm glass capillary ($\theta \approx 0$):

$$h = \frac{2\times0.072}{998\times9.81\times50\times10^{-6}} = 0.29\ \text{m}$$

**29 centimetres**, in a channel you can't see without a microscope. This is a
good demonstration of just how strong capillarity is at this scale, and it's
`code/ch10/capillary_rise.m`.

The crossover length is the **capillary length**:

$$\ell_c = \sqrt{\frac{\sigma}{\rho g}} \approx 2.7\ \text{mm for water}$$

Below $\ell_c$, surface tension dominates gravity — which for microfluidics means
**always**. $Bo = (L/\ell_c)^2 \approx 10^{-3}$ at 100 µm. Drop the gravity term.

---

## 10.5 The numbers that classify interface problems

### Capillary number — the important one

$$Ca = \frac{\mu U}{\sigma}$$

viscous shear vs surface tension. Use the **continuous-phase** viscosity — the
phase doing the shearing.

| $Ca$ | Regime | Interface behaviour |
|---|---|---|
| $<10^{-3}$ | capillary-dominated | Spheres/plugs. Squeezing in a T-junction. |
| $10^{-3}$–$10^{-2}$ | transitional | Dripping. |
| $>10^{-2}$ | viscous-dominated | Jetting, tip-streaming, threads. |

Tutorial 12 maps droplet-generation regimes onto exactly this axis.

### Weber and Ohnesorge

$$We = \frac{\rho U^2 L}{\sigma} = Re\cdot Ca, \qquad Oh = \frac{\mu}{\sqrt{\rho\sigma L}} = \frac{\sqrt{We}}{Re}$$

$We \ll 1$ at the microscale (droplets don't splash). $Oh$ compares viscous
damping to inertia–capillary oscillation; $Oh > 1$ means a pinching thread is
overdamped and breaks slowly rather than snapping.

---

## 10.6 Computing curvature — the thing that actually matters numerically

Tutorials 11 and 12 need $\kappa$ on a grid, from an implicit interface function
$\phi$ (level set, §11.2). The formula:

$$\kappa = \nabla\cdot\left(\frac{\nabla\phi}{|\nabla\phi|}\right)$$

In 2D, expanded:

$$\kappa = \frac{\phi_{xx}\phi_y^2 - 2\phi_x\phi_y\phi_{xy} + \phi_{yy}\phi_x^2}{(\phi_x^2+\phi_y^2)^{3/2}}$$

```matlab
function k = curvature(phi, dx, dy)
[px, py]   = grad2(phi, dx, dy);
[pxx, pxy] = grad2(px,  dx, dy);
[~,   pyy] = grad2(py,  dx, dy);
den = (px.^2 + py.^2).^1.5;
k   = (pxx.*py.^2 - 2*px.*py.*pxy + pyy.*px.^2) ./ max(den, eps);
end
```

**Curvature is a second derivative, so it amplifies noise by $1/\Delta x^2$.**
Three practical consequences, all of which bite in Tutorial 11:

1. Keep $\phi$ smooth — this is *why* level sets are reinitialized to a signed
   distance function (§11.3). A noisy $\phi$ gives a noisy $\kappa$ gives a noisy
   force gives spurious currents.
2. Curvature error drives **spurious (parasitic) currents** — unphysical
   vortices at the interface that never go away, because the discrete surface
   tension force is not exactly balanced by the discrete pressure gradient.
3. Resolution: you need **≥ 10 cells per droplet radius** for curvature to be
   reasonable, and 20+ if you care about the pressure jump.

`code/ch10/young_laplace.m` verifies the curvature routine against a circle
(where $\kappa = 1/R$ exactly) and reports the convergence rate — do this before
trusting any two-phase result.

---

## 10.7 Static meniscus shapes

For a 2D meniscus in a channel of width $w$ with contact angle $\theta$ on both
walls, the equilibrium interface is a **circular arc** (constant curvature, since
$\Delta p$ is constant and gravity is negligible at $Bo\ll1$), meeting each wall
at $\theta$. Geometry gives

$$R = \frac{w}{2\cos\theta}, \qquad \Delta p = \frac{\sigma}{R} = \frac{2\sigma\cos\theta}{w}$$

That $\Delta p$ is the pressure needed to push an interface through a
constriction of width $w$ — the **capillary entry pressure**, and it is the
number that governs:

- Whether liquid spontaneously wicks into your channel ($\theta<90°$, it pulls
  itself in) or must be forced ($\theta>90°$).
- Burst-valve design: a sudden expansion creates a curvature barrier that stops
  flow until a threshold pressure is applied. Passive valving in paper and
  centrifugal microfluidics works exactly this way.
- Why a narrow constriction blocks a non-wetting phase.

`code/ch10/meniscus.m` computes arc shapes and entry pressures across $\theta$
and $w$, and includes a burst-valve calculation.

---

## Exercises

**10.1** Run `young_laplace`. Confirm the computed curvature of a circle
converges to $1/R$ at second order. Then add random noise of amplitude
$10^{-3}$ to $\phi$ and recompute. How much does $\kappa$ degrade? *(This is
§10.6's warning made concrete, and it motivates reinitialization.)*

**10.2** Compute $\Delta p$ for water droplets of 1 µm, 10 µm, 100 µm, 1 mm
diameter. Compare with the ~0.7 Pa viscous pressure drop of Tutorial 8's channel.
At what droplet size do they become comparable?

**10.3** Run `capillary_rise`. Verify the $1/r$ scaling. At what radius does the
rise height equal the capillary length?

**10.4** Run `meniscus`. Compute the entry pressure for water ($\theta = 105°$,
native PDMS) into a 20 µm constriction. Is this achievable with a standard
syringe pump?

**10.5 (burst valve)** Design a passive burst valve: a hydrophilic channel that
wicks spontaneously, then a sudden expansion that halts the flow until 50 mbar is
applied. What expansion ratio do you need?

**10.6** For a 50 µm water droplet in oil ($\mu_{oil} = 10$ mPa·s,
$\sigma = 5$ mN/m) at 1 mm/s: compute $Ca$, $We$, $Oh$, $Bo$. Predict the regime
from the §10.5 table. You will test this prediction in Tutorial 12.

**10.7** Derive the capillary rise formula yourself from a force balance on the
liquid column ($2\pi r\sigma\cos\theta$ upward vs $\rho g \pi r^2 h$ downward).
Then explain why the 2D channel version has $w$ instead of $2/r$ and gives a
different numerical coefficient.

---

**Next:** [Tutorial 11 — Two-phase flow with level sets](11-level-set-two-phase.md).
Now you make the interface *move*.
