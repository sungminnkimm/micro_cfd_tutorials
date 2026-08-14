# Tutorial 14 — From a Question to a Model

**Goal:** The meta-skill. You have the tools; this is how to point them at a
problem nobody handed you.

---

## 14.1 The workflow

### Step 1 — State the question as a number

Not "how does my mixer perform" but **"what channel length gives $M<0.05$ for
fluorescein at 2 mm/s in a 50 µm channel?"**

A CFD project without a target number produces pretty pictures and no decisions.
Write the question down, with units, before you open MATLAB.

### Step 2 — Scaling analysis (Tutorial 1)

```matlab
s = scaling('L',50e-6,'U',2e-3,'D',4.3e-10);
```

This is free and it usually settles the modelling decisions:

| If | Then |
|---|---|
| $Re < 1$ | Drop inertia. Flow is steady → solve once, freeze it |
| $Pe > 10$ | Advection-dominated scalar → TVD scheme mandatory |
| $Bo < 0.01$ | Drop gravity |
| $Ca < 0.01$ | Interfaces are capillary-dominated → squeezing regime |
| $Kn > 0.001$ | Continuum assumption failing (gases below ~1 µm only) |

### Step 3 — Choose the simplest model that can answer the question

**Always ask whether you need CFD at all.** In rough order of cost:

1. **Algebraic.** Poiseuille, cubic law, $L_{mix}\approx0.3\,Pe\,w$,
   Young–Laplace, Garstecki. Answers a surprising number of design questions in
   seconds.
2. **1D transient.** The Tutorial 4.5 trick: a steady 2D problem at high $Pe$
   becomes an unsteady 1D one via $x = Ut$. Exact when axial diffusion is
   negligible.
3. **2D CFD.** Structure, recirculation, geometry effects.
4. **3D CFD.** When the answer genuinely depends on the third dimension —
   herringbone mixers, Dean flow, corner leakage in droplet squeezing.

Escalate only when the level below cannot answer the question. Tutorial 4's
Y-junction study is the case in point: the 1D model gave the mixing length
directly and revealed that the $Pe\,w$ rule overestimates by 2.9×.

### Step 4 — Resolution budget (Tutorials 5, 8)

Three constraints, all mandatory:

```matlab
dx_geom = smallest_feature/20;
dx_grad = channel_width/20;
dx_num  = 0.2*D/U;              % Pe_cell < 0.2 for upwind
```

With a TVD limiter, `dx_num` relaxes by ~50–100× (Tutorial 5 measured 757× less
$D_{num}$ for van Leer on a smooth profile). Then estimate the cost:

```matlab
dt = 0.4*min([dx/U, 0.5/nu/(2/dx^2), sqrt(rho*dx^3/(4*pi*sigma))]);
nsteps = tEnd/dt;
```

**If `nsteps > 1e6`, stop and rethink** before spending the CPU. `regime_map.m`
does this arithmetic for droplet problems specifically.

### Step 5 — Build up, never all at once

1. Straight channel → verify Poiseuille, `divMax`, mass balance.
2. Add geometry → re-verify.
3. Add scalar transport → measure $D_{num}/D$.
4. Add two phases → run `static_droplet` at your properties first.

Each step gets its own verification. Tutorial 13.7 lists three real bugs from
this course, every one of which produced plausible output and was caught only by
a diagnostic with a known right answer.

### Step 6 — Grid convergence, always

Two grids minimum, three for a GCI. If the answer moves, you haven't converged.

### Step 7 — Write down what you did *not* model

The professional habit. 2D? Static contact angle? Constant properties? No
surfactant dynamics? Say so.

---

## 14.2 Non-dimensionalization: worth the effort

Running in SI at the microscale means numbers like $10^{-6}$ m and
$10^{-10}$ m²/s, which is fine numerically but makes errors hard to spot.
Non-dimensionalizing gives $O(1)$ numbers where a wrong answer *looks* wrong.

Scale by $L$, $U$, $L/U$, and $\mu U/L$ (the **viscous** pressure scale — see
Tutorial 1.2; using $\rho U^2$ at $Re\ll1$ is wrong by a factor of $1/Re$):

$$Re\left(\frac{\partial\mathbf{u}^*}{\partial t^*} + \mathbf{u}^*\cdot\nabla^*\mathbf{u}^*\right) = -\nabla^*p^* + \nabla^{*2}\mathbf{u}^*$$

Now one simulation at $Re = 0.1$ covers **every** device with that $Re$.

The catch: with two phases you need $Ca$ *and* $Re$ *and* density and viscosity
ratios, and matching all of them is harder than just running in SI. This
course's solvers are dimensional for that reason, and because unit errors are
easier to catch when the numbers carry units you recognize.

---

## 14.3 Extending the code

### Electrokinetics (electroosmotic flow)

Not covered in the course; here is how to bolt it on. In the thin-EDL limit
(Debye length $\lambda_D\sim$ nm $\ll$ channel), you do **not** resolve the
double layer. Instead apply the **Helmholtz–Smoluchowski slip velocity** at the
wall:

$$u_{slip} = -\frac{\varepsilon\zeta}{\mu}E_t$$

with $\zeta$ the zeta potential and $E_t$ the tangential field.

1. Solve $\nabla^2\varphi = 0$ for the applied potential (your `poisson2d`,
   Dirichlet at electrodes, Neumann at walls).
2. Compute $\mathbf{E} = -\nabla\varphi$.
3. Set the wall BC to `'moving'` with $u = u_{slip}$ — the solver already
   supports it.

The striking result: in a uniform channel this gives a **plug-flow** profile,
not parabolic. That's why capillary electrophoresis has so much less band
broadening than pressure-driven flow — no shear means no Taylor dispersion
(Tutorial 9.5). Predicting that from the two numbers is a good check that you've
understood both tutorials.

### Slip flow (gases below ~1 µm)

For $10^{-3} < Kn < 10^{-1}$, replace no-slip with Maxwell slip:

$$u_{wall} = \frac{2-\sigma_v}{\sigma_v}\lambda\frac{\partial u}{\partial n}$$

In `tang_ghost`, blend between the no-slip ghost ($-u_{in}$) and the free-slip
ghost ($+u_{in}$) according to $Kn$. A ten-line change.

### Non-Newtonian fluids

Blood, polymer solutions, many biological samples. Replace the constant `nu`
with $\mu(\dot\gamma)$ evaluated from the local strain rate:

$$\dot\gamma = \sqrt{2\left(\frac{\partial u}{\partial x}\right)^2 + 2\left(\frac{\partial v}{\partial y}\right)^2 + \left(\frac{\partial u}{\partial y}+\frac{\partial v}{\partial x}\right)^2}$$

with, e.g., Carreau–Yasuda. `ns2phase.m` already carries variable $\mu$ on the
faces, so that machinery is in place.

### Particles and cells

Point particles: advect with the local velocity plus a Stokes drag correction.
Finite-size particles: needs a proper immersed-boundary method — beyond a
straightforward extension.

### 3D

The numerics are identical; the cost is not. $256^3 = 1.7\times10^7$ unknowns is
beyond a direct Poisson factorization, so you'd need a multigrid or FFT-based
solver. Consider whether the question really requires it (Tutorial 6.5 has the
2D-error table) before committing.

---

## 14.4 A debugging playbook

| Symptom | Look here first |
|---|---|
| `NaN` in one step | Division by zero — an unguarded limiter ratio (Tutorial 5.5), or `1/\|∇φ\|` |
| Blows up over ~10 steps | Timestep. Check *all four* limits, not just CFL |
| Cell-to-cell sawtooth, growing | Stability violation (Tutorial 4.2) |
| Cell-to-cell wiggles, steady | $Pe_{cell}>2$ with central differencing (Tutorial 5.4) |
| Smooth but wrong | Boundary condition. Almost always |
| Converges at order 1, not 2 | Boundary condition (Tutorial 2.3) |
| `divMax` not at round-off | Pressure solve inconsistent with the corrector (Tutorial 8.1) |
| Mass not conserved | Outflow BC, or level-set reinit drift (Tutorial 11.3) |
| Result changes with grid | Not converged. This is not a bug, it's a measurement |
| Mixing looks too good | Numerical diffusion (Tutorial 5.6) — measure it |
| Interface wobbles | Spurious currents (Tutorial 11.4) — measure them |
| Droplet dissolves | Level-set mass loss — turn on the area correction |

**The meta-rule:** find the *smallest* case that reproduces the problem. A bug
in a 400×100 droplet simulation is unfindable; the same bug in a 8×8 static case
is obvious. Every solver in this course prints diagnostics for exactly this
reason.

---

## 14.5 Project templates

Each is sized to be genuinely completable, and each ends with a claim you could
defend.

**A. Mixer design.** Given two streams and a length budget, compare straight,
serpentine, and a lamination inlet. Deliverable: mixing length for each, with
measured $D_{num}/D$ and a GCI. *(Tutorials 8, 9, 13.)*

**B. Droplet generator.** Design for a target size and rate. Deliverable: the
regime map placement, a simulated size, comparison with Garstecki, and a grid
convergence study. *(Tutorials 10, 11, 12.)*

**C. Dispersion in a separation channel.** How much does Taylor–Aris broaden
your band? Compare pressure-driven with electroosmotic (using §14.3). *(Tutorial 9.)*

**D. Heat management.** Conjugate heat transfer in a microchannel heat sink:
what flow rate keeps the chip below 85 °C? *(Tutorials 3, 8, plus a temperature
scalar.)*

**E. Reproduce a paper.** Pick a microfluidics paper with a figure you can
extract numbers from. Reproduce it. Deliverable: an overlay plot and an honest
account of every discrepancy. *(This is the most valuable one.)*

---

## 14.6 What you should now be able to do

- Decide, from scaling, which terms matter and which model to use
- Discretize a PDE and **prove** the discretization is right
- Write a Navier–Stokes solver from scratch on a staggered grid
- Build your own geometry with inlets, outlets, and obstacles
- Transport scalars without lying to yourself about numerical diffusion
- Track interfaces and apply surface tension
- Predict droplet size and mixing length
- Quantify your own numerical error, and state what your model omits

That last one is what separates CFD from colouring in. Every tutorial in this
course printed a number with a known right answer — and three times, that number
caught a real bug that produced perfectly plausible output. Keep the habit.

---

## Appendices

- [A. MATLAB performance for CFD](A-matlab-performance.md)
- [B. Common bugs and their symptoms](B-common-bugs.md)
- [C. Notation and symbols](C-notation.md)
