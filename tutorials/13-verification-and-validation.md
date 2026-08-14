# Tutorial 13 — Verification and Validation

**Goal:** Prove your code is right. Not "it looks plausible" — prove it, with
numbers a reviewer can check.

**Code:** `code/ch13/mms_driver.m`, `code/ch13/gci.m`

---

## 13.1 The two questions

Roache's distinction, and it is worth internalizing because people conflate them
constantly:

> **Verification:** *Am I solving the equations right?*
> A pure mathematics question. Compare against exact solutions of the discrete
> or continuous problem. No physics involved.
>
> **Validation:** *Am I solving the right equations?*
> A physics question. Compare against experiment.

A code can be perfectly verified and completely invalid — a beautifully
converged 2D simulation of a device whose behaviour is 3D. Tutorial 6.5 measured
exactly that: the 2D flow-rate prediction for a square channel is 137% too high,
and no amount of grid refinement fixes it, because the error is in the *model*,
not the *discretization*.

**Verify first.** You cannot learn anything from comparing a buggy code to
experiment — the errors are indistinguishable.

---

## 13.2 The verification hierarchy

Work up. Each level assumes the ones below it pass.

| Level | What it checks | This course |
|---|---|---|
| 0 | Does it run? | — |
| 1 | Conservation (mass, momentum) | `divMax`, `massErr`, area error |
| 2 | Symmetry / invariance | fore-aft symmetry at $Re\ll1$ |
| 3 | Analytic solutions | Poiseuille, capillary rise, Taylor–Aris |
| 4 | **Order of accuracy (MMS)** | `mms_driver.m` |
| 5 | Benchmark comparison | Ghia cavity |
| 6 | Code-to-code comparison | — |

**Level 4 is the strongest.** Levels 1–3 can pass with a bug that happens to
preserve the tested property. Order-of-accuracy verification will not: a
first-order bug anywhere in a second-order code drags the *whole* observed rate
to 1, and it does so reliably.

### Measured results from this course

Every claim below is output from the code in this repository:

| Test | Result |
|---|---|
| FD stencils (T2) | orders 1.00 / 2.00 / 2.00 |
| Poisson MMS (T3) | order 1.99 |
| Time integration (T4) | FE 1.00, BE 1.00, CN 2.00 |
| Upwind $D_{num}$ (T5) | measured 6.250e-04 vs theory 6.250e-04 |
| Poiseuille (T6) | order 2.00; $u_{max}/\bar u = 1.498$ |
| Cavity vs Ghia, $Re=100$ (T7) | $\|u-u_{Ghia}\| < 0.005$; `divMax` 7e-15 |
| Channel (T8) | `divMax` 4e-14; $\Delta p$ within 0.09% |
| Taylor–Aris (T9) | $D_{eff}$ within 1.4% of theory |

That table is what "verified" looks like. Assemble one for your own code.

---

## 13.3 The Method of Manufactured Solutions

The most powerful verification tool, and it works for *any* PDE — including ones
with no known exact solutions.

**The trick:** stop looking for solutions. Pick one, then work out what source
term would make it true.

1. **Choose** a solution $\mathbf{u}_{ex}, p_{ex}$. It need not be physical —
   just smooth and non-trivial. Trigonometric functions are ideal.
2. **Substitute** into your PDE. It won't balance; the residual is a source
   term $\mathbf{f}$:
   $$\mathbf{f} = \frac{\partial\mathbf{u}_{ex}}{\partial t} + (\mathbf{u}_{ex}\cdot\nabla)\mathbf{u}_{ex} + \frac{1}{\rho}\nabla p_{ex} - \nu\nabla^2\mathbf{u}_{ex}$$
3. **Add $\mathbf{f}$** to your solver as a body force, with BCs from $\mathbf{u}_{ex}$.
4. **Solve.** Your code must now reproduce $\mathbf{u}_{ex}$.
5. **Refine** and measure the order.

Do the differentiation with the Symbolic Toolbox if you have it, or by hand.
`mms_driver.m` does it by hand for a divergence-free trigonometric field so no
toolbox is needed.

### Choosing the manufactured solution

- **Make it divergence-free** if your solver assumes incompressibility. Use a
  streamfunction: $u = \partial\psi/\partial y$, $v = -\partial\psi/\partial x$
  is divergence-free automatically.
- **Avoid symmetries** that accidentally cancel error terms. $\sin(2\pi x)$ on
  $[0,1]$ has a symmetry that can hide bugs; $\sin(2\pi x + 0.3)$ doesn't.
- **Exercise every term.** If $\mathbf{u}_{ex}$ is linear in $x$, the convective
  term is trivially satisfied and you have verified nothing about it.
- **Non-square grid.** $n_x \ne n_y$ catches transposition bugs instantly
  (Tutorial 3.5) — this course's Poisson MMS deliberately uses $n_y = 3n_x/4$.

### Reading the result

| Observed order | Diagnosis |
|---|---|
| 2.0 | Correct |
| 1.0 | First-order bug — usually a boundary condition |
| 1.5 | Mixed: one boundary or one term is first order |
| 0 | Not converging — a real error, not a discretization error |
| 3+ | Suspicious. Your manufactured solution is probably too simple. |

---

## 13.4 Grid convergence for problems with no exact answer

Once you're solving a real problem, there's no exact solution to compare to. Use
**Richardson extrapolation** to estimate the exact answer from three grids.

With refinement ratio $r$ (usually 2) and solutions $f_1$ (finest), $f_2$, $f_3$:

$$p = \frac{\ln\left|\frac{f_3-f_2}{f_2-f_1}\right|}{\ln r}, \qquad f_{exact} \approx f_1 + \frac{f_1-f_2}{r^p-1}$$

The **Grid Convergence Index** turns this into an error bar:

$$GCI = \frac{F_s\,|f_1-f_2|}{(r^p-1)\,|f_1|}\times100\%$$

with $F_s = 1.25$ for three grids (a safety factor). GCI is the standard way to
report numerical uncertainty in a paper — `gci.m` computes it.

**Report it.** "The predicted mixing length was 3.4 mm with a grid convergence
index of 2.1%" is a defensible claim. "We used a fine mesh" is not.

### The check people skip

$p$ from the formula above should come out near your *theoretical* order. If you
get $p = 0.4$ or $p = 6$, the grids are not in the **asymptotic range** and
Richardson extrapolation is invalid — the error bar you'd compute is
meaningless. Refine further before believing it.

---

## 13.5 Validation, and the honest report

Once verified, compare with experiment. Sources of discrepancy, in the order you
should suspect them:

1. **Model form.** 2D vs 3D (the big one — Tutorial 6.5), missing physics,
   wrong constitutive law.
2. **Input uncertainty.** Do you actually know $\sigma$ for your surfactant
   system? Channel depth to ±5%? Fabrication tolerances are usually the
   dominant uncertainty in microfluidics — and recall the cubic law, where a
   10% depth error is a 30% flow-rate error.
3. **Boundary conditions.** Is your inlet really a parabola? Is the outlet
   really at zero gauge pressure?
4. **Numerical error.** You already quantified this. It's usually the smallest.

That ordering is the point. Beginners assume discrepancy means numerical error,
and grind on grid refinement while a 137% model error sits untouched.

### What a defensible claim looks like

> "The predicted droplet length was $L/w = 1.48$, compared with $1.52\pm0.03$
> measured. The grid convergence index was 1.8% (three grids, observed order
> 1.9). Level-set mass loss over the run was 0.4%. The simulation is 2D and
> therefore does not capture corner leakage of the continuous phase, which is
> expected to reduce $L/w$ by a further few percent."

Every number has a provenance. Nothing is asserted that wasn't measured. The
known model limitation is stated rather than hoped over.

---

## 13.6 A verification checklist for your own code

Copy this into your project:

```
[ ] Conservation
    [ ] max|div u| ~ 1e-14 after projection
    [ ] mass balance in = out to round-off
    [ ] scalar mass conserved (closed domain)
    [ ] level-set area drift measured and reported
[ ] Order of accuracy
    [ ] spatial order measured by MMS on a NON-SQUARE grid
    [ ] temporal order measured with the grid held fixed
    [ ] both match theory
[ ] Analytic benchmarks
    [ ] Poiseuille profile and pressure drop
    [ ] at least one problem exercising every term you use
[ ] Physical sanity
    [ ] symmetry where physics demands it
    [ ] correct limiting behaviour (Re -> 0, Ca -> 0)
[ ] Numerical error budget
    [ ] D_num/D_phys measured for scalar transport
    [ ] spurious currents vs physical velocity for two-phase
    [ ] GCI reported for the quantity of interest
[ ] Model limitations stated
    [ ] 2D vs 3D
    [ ] missing physics named explicitly
```

---

## 13.7 A worked case: the bugs this course actually hit

Verification is not hypothetical. Three real failures were caught while building
this material, each by a different level of the hierarchy — and each looked
plausible until measured.

**1. Masked pressure operator (Level 1, conservation).** The contraction
geometry ran, produced a sensible-looking velocity field, and gave the right
throat velocity (2.001× upstream, exactly as mass conservation demands). But
`divMax` was **416**. The mask had been applied to velocities and not to the
Poisson operator, so the projection was making the flow divergence-free with
respect to a domain that had no obstacle in it. *Caught only because `divMax`
was printed every run.*

**2. Outflow face not pressure-corrected (Level 1).** After fixing (1), `divMax`
was 1.5 — better, but not round-off. The outlet velocity was being extrapolated
after the projection instead of corrected with the Dirichlet-$\phi$ ghost the
operator had assumed. The error lived entirely in the last column of cells.
*Caught because "1.5 is not 1e-14" and the distinction was treated as
meaningful.*

**3. Curvature sign (Level 3, analytic).** The static droplet gave a pressure
jump of **−169 Pa where theory says +125 Pa**. With $\phi>0$ inside, $\nabla\phi$
points inward, so the standard formula returns $\kappa=-1/R$, and surface
tension was pushing the droplet apart instead of holding it together. *Caught
because there is an exact answer, $\sigma/R$, and it was checked rather than
eyeballed.* Fixing it dropped spurious currents by ~30 orders of magnitude.

**4. Pressure drift at a outlet (Level 1) — plus a wrong diagnosis.** The
serpentine diverged at step 167 with a timestep at *half* the viscous limit,
while `max|div u|` sat at $10^{-13}$ throughout. The first hypothesis was a
cell-Reynolds instability, from a mental estimate of $Re_{cell}\approx7.5$.
Dropping the convective term made it survive 36× longer, which felt like
confirmation — and was not. The estimate was off by 1000×; the true
$Re_{cell}$ was 0.0075.

Printing *where* the blow-up lived settled it in one line: the outlet face, with
$p\sim10^{154}$. In the incremental scheme $p^{n+1}=p^n+\rho\phi$, so setting
$\phi=0$ at a pressure outlet pins the *correction*, not the *pressure*, and the
accumulated pressure drifts without bound. The fix is $\phi_{wall}=-p_{wall}/\rho$.

### What these have in common

**All four produced output that looked reasonable.** None crashed initially. Two
gave a *correct* answer to the first quantity anyone would check. They were
caught by diagnostics that print a number with a known right answer — which is
the entire argument for building those diagnostics in from the start.

Case 4 adds a second lesson, about the debugging itself:

- **A healthy `divMax` certifies the projection and nothing else.** It says
  nothing about whether the momentum equation or the pressure is sane.
- **Locate before you theorize.** `[~,k] = max(abs(u(:)))` followed by
  `ind2sub` took seconds and pointed straight at the outlet. An hour went into a
  plausible theory that one `fprintf` of the actual number would have killed.
  When your hypothesis contains a number, print the number.
- **A partial improvement is not confirmation.** Dropping convection helped
  measurably and was still the wrong fix. "It got better" is weak evidence; "it
  is now at round-off" is strong evidence.

---

## Exercises

**13.1** Run `mms_driver`. Confirm second order for the full NS solver.

**13.2** Introduce a deliberate bug: change one ghost-cell formula from
`2*u_wall - u_in` to `-u_in`. Rerun the MMS. What order do you measure?
*(This is the single most valuable exercise in the tutorial — see a bug caught
by a rate rather than by inspection.)*

**13.3** Run `gci` on the mixing length from Tutorial 9, using three grids.
Report the GCI. Are you in the asymptotic range?

**13.4** Take one result from Tutorial 12 and write the §13.5 paragraph for it.
Every number must have a provenance.

**13.5** Build the §13.6 checklist for your own project and work through it. Note
which boxes you cannot yet tick — that list is your next week of work.

---

**Next:** [Tutorial 14 — From a question to a model](14-your-own-model.md).
