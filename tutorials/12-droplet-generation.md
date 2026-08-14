# Tutorial 12 — Droplet Generation

**Goal:** The capstone. Simulate a T-junction and a flow-focusing device,
identify the operating regime, and predict droplet size. This uses every piece
of the course at once.

**Code:** `code/ch12/t_junction.m`, `code/ch12/flow_focusing.m`, `code/ch12/regime_map.m`

---

## 12.1 Why droplets

Each droplet is an independent picolitre reactor. Generate them at kHz rates and
you get:

- **Digital PCR** — one DNA molecule per droplet, count the positives
- **Single-cell sequencing** — one cell per droplet with a barcoded bead
- **Directed evolution** — one variant per droplet, screened by fluorescence
- **Crystallization / kinetics screening** — thousands of conditions per chip

The requirement is always the same: **monodispersity**. Coefficient of variation
below ~2%, or the statistics that make the assay work fall apart. Which makes
"what sets the droplet size, and how sensitive is it?" the central design
question — and it's what this tutorial computes.

---

## 12.2 The T-junction

Dispersed phase enters through a side channel; continuous phase flows down the
main channel. Three regimes, set mainly by the **capillary number**
$Ca = \mu_c U_c/\sigma$ (continuous-phase viscosity — the phase doing the
shearing).

### Squeezing ($Ca \lesssim 10^{-2}$)

The dominant regime in practice, and it is not what intuition suggests.

The dispersed phase grows until it **blocks the main channel**. Blocking forces
the continuous phase through a thin gap, and the resulting *pressure* build-up —
not shear — squeezes the neck until it pinches off.

Because the mechanism is pressure/geometry rather than shear, the droplet size
is essentially **independent of $Ca$**, and given by the Garstecki scaling law
(Garstecki et al., *Lab Chip* 2006):

$$\boxed{\frac{L}{w} = 1 + \alpha\frac{Q_d}{Q_c}}$$

with $\alpha \approx 1$ (geometry-dependent, typically 1–2). $L$ is the plug
length, $w$ the main-channel width, $Q_d/Q_c$ the flow-rate ratio.

**This is the most useful formula in droplet microfluidics.** It says droplet
size is set by the *flow-rate ratio* and *channel geometry* alone — not by
viscosity, not by surface tension, not by absolute flow rate. That's why the
squeezing regime is so reproducible and why devices are designed to sit in it.

### Dripping ($10^{-2} \lesssim Ca \lesssim 10^{-1}$)

Shear becomes comparable to interfacial tension. Droplets pinch off before
fully blocking the channel, and are **smaller than the channel width**. Size now
depends on $Ca$:

$$\frac{d}{w} \sim Ca^{-1/3} \quad\text{(approximately)}$$

More tunable, but also more sensitive to fluctuations in viscosity, temperature,
and surfactant — hence usually worse monodispersity.

### Jetting ($Ca \gtrsim 10^{-1}$)

The dispersed phase forms a long jet that breaks up downstream via the
**Rayleigh–Plateau** instability. Poor monodispersity. Usually avoided — but note
that Rayleigh–Plateau is fundamentally 3D, so a 2D simulation will not reproduce
jetting correctly.

---

## 12.3 Flow focusing

Dispersed phase in the centre, continuous phase from both sides, everything
forced through a narrow orifice. Advantages: smaller droplets (below the channel
width), higher rates, better control.

Scaling in the dripping regime:

$$\frac{d}{w_{or}} \sim \left(\frac{Q_d}{Q_c}\right)^{1/3}Ca^{-1/3}$$

The orifice width is the controlling dimension, which is why flow-focusing
devices can make droplets far smaller than their channels.

---

## 12.4 What you must get right numerically

Droplet generation is the most demanding simulation in this course. Checklist:

**1. Resolution.** ≥ 20 cells across the neck at pinch-off, and the neck gets
*thin*. This is the binding constraint and it is expensive. If you under-resolve
the neck, breakup happens when the level set decides two interfaces are close
enough — which is a **numerical** criterion, not a physical one, and your droplet
size will be a function of $\Delta x$.

**2. The capillary timestep.** From Tutorial 11:
$\Delta t \le \sqrt{\bar\rho\Delta x^3/(4\pi\sigma)}$. This dominates. Budget for
$10^5$–$10^6$ steps.

**3. Spurious currents below the physical velocity.** Check with
`static_droplet` at your grid and fluid properties *before* running a
generation case. If parasitic velocities are comparable to $U_c$, your droplet
sizes are noise.

**4. Mass conservation.** Level sets leak. Report the area error. If you lose
5% of your droplet volume, your size prediction is 5% wrong — which for a
technique whose selling point is 2% CV is fatal.

**5. Wetting.** The continuous phase must wet the walls. This code doesn't model
contact lines (Tutorial 11.7), which is *fine* for a properly wetted device
because the dispersed phase never touches a wall — but it means the code cannot
predict what happens when wetting fails, which is the most common real-world
failure mode.

### The honest 2D caveat

A 2D "droplet" is an infinite cylinder. So:

- Volume scaling with size differs ($L^2$ vs $L^3$).
- The pressure that drives squeezing differs, because in 3D the continuous
  phase leaks through the **corners** of a rectangular channel — a genuinely 3D
  effect that is central to the squeezing mechanism.
- Rayleigh–Plateau breakup is absent.

**2D gets the regimes and the trends right; it gets absolute droplet volumes
wrong.** Use it to understand *how* your device behaves and which direction to
move a parameter. Use experiments, or 3D, for absolute numbers.

---

## 12.5 The recipe

```matlab
regime_map                % where does my design sit? (do this FIRST, it's free)
t_junction                % simulate it
flow_focusing             % compare geometries
```

`regime_map.m` needs no simulation at all — it evaluates $Ca$, $We$, $Bo$, the
Garstecki prediction, and the required timestep and step count for a given grid.
**Run it before committing hours of CPU.**

### A cautionary tale, from building this course

The first draft of `t_junction.m` used $L_x = 500$ µm on a 200×80 grid with a
10 mPa·s oil at 2 mm/s. Perfectly reasonable-looking physical parameters. The
arithmetic:

$$\Delta x = 2.5\ \mu\text{m} \Rightarrow \Delta t_{visc} = 1.4\times10^{-7}\ \text{s},\quad \Delta t = 5.6\times10^{-8}\ \text{s}$$
$$t_{res} = \frac{500\ \mu\text{m}}{2\ \text{mm/s}} = 0.3\ \text{s} \Rightarrow \boxed{5,333,333\ \text{steps}}$$

about 90 minutes for one case. It was discovered by watching a test suite fail
to finish, not by calculating — even though `regime_map.m`, written for exactly
this purpose two files earlier, would have reported it in one second.

The fix illustrates the trade-offs worth knowing:

| Change | Effect on $\Delta t$ | Effect on $Ca$ |
|---|---|---|
| Coarsen the grid | $\Delta t_{visc}\propto\Delta x^2$, $\Delta t_{cap}\propto\Delta x^{3/2}$ — **big win** | none |
| **Lower $\mu_c$** | $\nu$ falls → $\Delta t_{visc}$ rises — **win** | $Ca$ falls — **also a win** |
| Raise $U_c$ | shortens $t_{end}$ — win | $Ca$ rises — pushes toward dripping |
| Lower $\sigma$ | $\Delta t_{cap}$ rises — win | $Ca$ rises — toward dripping |

Lowering the continuous-phase viscosity helps *twice*, which is why the shipped
version uses 2 mPa·s. The shipped configuration runs in ~40k steps — but at only
~17 cells across the channel, **below** the ≥20 that §12.4 demands. That is
stated in the file's header rather than hidden, and it is why Exercise 12.5 (grid
convergence) is not optional.

`ns2phase.m` now refuses to start a run exceeding `opt.maxSteps` (default
$3\times10^5$) and prints these options, so the mistake fails loudly instead of
quietly consuming an afternoon.

`t_junction.m` runs the simulation and measures droplet length, comparison with
Garstecki, and the diagnostics (area error, `divMax`, cells across the channel).

### What the shipped T-junction actually does — read this before running it

**It does not produce detached droplets.** That is the honest measured outcome at
the affordable configuration, and it is worth more to you than a prettier one.

The run is numerically healthy — `max|div u| ≈ 2×10⁻¹²`, and the enclosed
dispersed area *grows* (+0.47 relative) exactly as continuous injection demands.
The dispersed phase enters the junction, fills it, and is dragged downstream. But
at 16.7 cells across the channel, the **neck** is only a few cells wide, and §12.4
asks for ≥20 across the neck. Pinch-off is never resolved.

Now consider the alternative. If it *had* pinched off at this grid, breakup would
have been triggered by two interfaces coming within a cell of each other — a
**numerical** criterion, not a physical one — and the resulting droplet length
would have been a function of $\Delta x$ rather than of your device. **A demo that
silently produced a plausible droplet size at this resolution would be worse than
one that refuses to.** You would have believed it.

This is also why the area diagnostic matters. An earlier version of this script
inherited the default reinitialization settings (every 20 steps, 3 iterations).
Over 50,000 steps that is 2,500 reinit calls at ~0.1% area each — and the
interface *bled away faster than it was injected*, giving a **negative** area
error of −0.43 in a run where fluid is continuously added. Reinitializing every
60 steps instead flipped it to +0.47. A conservation diagnostic whose *sign* is
wrong is telling you something loudly; read it.

Exercise 12.5 is where you fix the resolution. Expect to need ~3–4× this grid,
and run `regime_map` for the step count before you commit.

---

## 12.6 Reading the results

| Symptom | Likely cause |
|---|---|
| Droplets never pinch off | Under-resolved neck, or $Ca$ too low for the domain length |
| Droplet size drifts over time | Mass loss in the level set |
| Size changes when you refine the grid | Not converged — the neck is under-resolved |
| Interface wobbles/oscillates | Spurious currents, or $\Delta t$ too near the capillary limit |
| Dispersed phase wets the wall | Missing contact-angle model; in a real device, a surface-treatment failure |
| Frequency but no clean droplets | Jetting regime — reduce $Ca$ |

The second and third rows are the ones to take seriously, because they produce
**plausible wrong answers** rather than obvious failures.

---

## 12.7 A worked design problem

*Make 50 µm water-in-oil droplets at 100 Hz, CV < 2%.*

1. **Choose the regime.** Monodispersity → squeezing → $Ca < 10^{-2}$.
2. **Choose fluids.** Oil $\mu_c = 10$ mPa·s, $\sigma = 5$ mN/m with surfactant.
   $Ca = \mu_c U_c/\sigma < 10^{-2}$ gives $U_c < 5$ mm/s.
3. **Choose geometry.** Squeezing gives plugs of length $\sim w$, so for 50 µm
   droplets use $w \approx 50$ µm.
4. **Flow rates.** $L/w = 1 + Q_d/Q_c$. For $L/w \approx 1.5$, $Q_d/Q_c = 0.5$.
5. **Check the frequency.** $f = Q_d/V_{drop}$. With $h = 50$ µm,
   $V \approx 50\times50\times75$ µm³ $\approx 190$ pL. For 100 Hz,
   $Q_d = 19$ nL/s = 1.1 µL/min, so $Q_c = 2.3$ µL/min.
6. **Verify $U_c$.** $Q_c/(wh) = 2.3\ \mu\text{L/min}/(50\times50\ \mu\text{m}^2)
   \approx 15$ mm/s → $Ca = 0.03$. **Too high** — that's dripping, not squeezing.
7. **Iterate.** Lower the flow rates and accept a lower frequency, or widen the
   channel, or use a less viscous oil.

**Step 6 is the point of the exercise.** The naive design is self-inconsistent,
and you only find out by checking. `regime_map.m` automates exactly this loop —
it takes a target and tells you whether the numbers close.

---

## Exercises

**12.1** Run `regime_map` for the §12.7 design. Confirm the inconsistency at
step 6. Find a self-consistent parameter set.

**12.2** Run `t_junction`. Measure the plug length and compare with
$L/w = 1 + \alpha Q_d/Q_c$. What $\alpha$ do you get?

**12.3** Vary $Q_d/Q_c$ over 0.25, 0.5, 1, 2. Confirm the linear Garstecki
relation. Plot $L/w$ vs $Q_d/Q_c$ and fit $\alpha$.

**12.4** Vary $Ca$ over $10^{-3}$, $10^{-2}$, $10^{-1}$ at fixed $Q_d/Q_c$.
Identify where the size stops being $Ca$-independent — that's your
squeezing→dripping transition. Compare with the §12.2 boundaries.

**12.5 (the convergence check that decides whether any of this is real)** Run
your best T-junction case on three grids. Plot droplet size vs $\Delta x$. If it
hasn't plateaued, you have measured your mesh, not your device. Report the
converged value *and* the grid you needed.

**12.6** Run `flow_focusing`. Compare droplet size to the T-junction at the same
$Ca$ and flow ratio. Which gives smaller droplets? Does it match §12.3?

**12.7** Add surfactant: drop $\sigma$ from 5 to 1 mN/m at fixed velocity. What
happens to $Ca$, to the regime, and to droplet size? *(This is the single most
common experimental knob — predict its effect before you turn it.)*

**12.8 (capstone)** Take your device from Exercise 1.5. Design a droplet
generator for it: pick geometry, fluids, and flow rates; predict size and
frequency analytically; simulate to check; run a grid-convergence study; and
write a one-page report stating what your 2D simulation does and does not
establish. That last part is the professional skill this course has been
building toward.

---

**Next:** [Tutorial 13 — Verification and validation](13-verification-and-validation.md).
You have built everything. Now learn to prove — to yourself and to a reviewer —
that it's right.
