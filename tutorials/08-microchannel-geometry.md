# Tutorial 8 — Building Your Own Geometry

**Goal:** Point the solver at a channel *you* designed. Inlets, outlets,
obstacles, contractions, serpentines. This is where the course stops being
exercises and starts being modelling.

**Code:** `code/ch08/channel.m`, `code/ch08/obstacle.m`, `code/ch08/serpentine.m`, `code/lib/make_mask.m`

---

## 8.1 Inflow and outflow: the compatibility trap

A cavity is closed. A channel is not, and open boundaries are where most
beginner solvers fail.

### Inlet

Prescribe the velocity profile directly. Because $u$ lives *on* the $x$-faces
(Tutorial 6), this is a direct assignment — no ghost cells:

```matlab
opt.bc.left = struct('type','inflow','u', uProfile);   % length ny
```

**Which profile?** Two choices, and the right one is usually the second:

1. **Uniform (plug) flow.** Physically corresponds to a reservoir. Requires an
   entrance length to develop.
2. **Parabolic (Poiseuille).** The fully-developed profile.

Use the parabola. Tutorial 1 computed the entrance length as $0.06\,Re\,L$, which
at $Re = 0.1$ is **0.6 µm** — utterly negligible. The flow is fully developed
essentially at the inlet, so imposing the parabola is not a shortcut, it's the
physically correct condition, and it saves you from wasting grid on a
development region that doesn't exist.

```matlab
yq = opt.Ly;   yc = G.yu;
uProfile = 6*Umean * (yc/yq) .* (1 - yc/yq);    % parabola, mean = Umean
```

Check: $\int_0^h 6\bar u \frac{y}{h}(1-\frac{y}{h})\,dy/h = \bar u$. ✓

### Outlet

Harder, because you don't know the answer there — that's what you're solving
for. Options, in increasing order of quality:

| Condition | Formula | Notes |
|---|---|---|
| Zero gradient | $\partial u/\partial x = 0$ | Simple, reflects a little |
| Convective | $\partial u/\partial t + U_c\partial u/\partial x = 0$ | Best for unsteady/vortices |
| Pressure outlet | $p = 0$, $\partial u/\partial n = 0$ | What we use |

We use the third: $\phi$ gets a **Dirichlet** condition (`'D'`) at the outlet in
the Poisson solve, which conveniently also removes the all-Neumann singularity
from Tutorial 3.

### The trap: correcting the outlet face

Here is the failure that bit *this* code, and the fix is not the obvious one.

Recall the compatibility condition (Tutorial 3, §3.3): the Poisson problem
$\nabla^2\phi = \nabla\cdot\mathbf{u}^*/\Delta t$ is solvable only if
$\oint \mathbf{u}^*\cdot\mathbf{n}\,dS = 0$. The tempting reading is "so I must
force inflow to equal outflow by rescaling the outlet profile."

**That reasoning applies only to an all-Neumann problem.** With a pressure
outlet, $\phi$ has a *Dirichlet* condition there, the operator is non-singular,
and no compatibility condition is required. Rescaling the outlet flux by hand is
solving a problem you don't have — and it papers over the one you do.

The real requirement is **consistency between the operator and the corrector.**
The Poisson operator assumed flux flows through the outlet face (that's where
the $-2c$ diagonal entry comes from, via the ghost $\phi_g = -\phi_c$). So the
corrector must apply the matching velocity correction on that face:

```matlab
% face sits half a cell from the centre, and phi_wall = 0
if bcT{2}=='D', u(end,:) = us(end,:) - dt*(0 - phi(end,:))*2/dx; end
```

If you instead re-extrapolate `u(end,:) = u(end-1,:)` after the projection, the
operator and the corrector disagree, and the last column of cells keeps a finite
divergence. **Measured in this code: `max|div u|` was 1.5 with re-extrapolation
and 4.2e-14 once the face was properly corrected.**

Global mass balance then comes out at round-off *as a consequence* of
$\nabla\cdot\mathbf{u}=0$ everywhere, rather than because it was imposed. That's
the right way round: a conservation property you had to enforce by hand is
usually a symptom that something upstream is inconsistent.

**Check it in your own runs:**

```matlab
fprintf('mass balance: %.3e   divMax: %.3e\n', S.massErr, S.divMax);
```

Both should be $\sim10^{-14}$ or better. If they aren't, stop and fix that
before believing anything else the solver tells you.

---

## 8.2 Validating the channel

Before adding any geometry, run a **plain straight channel** and check it
reproduces Poiseuille flow. `code/ch08/channel.m` does this and reports:

- Outlet profile vs the analytic parabola — should agree to the discretization
  error you measured in Tutorial 6
- $u_{max}/\bar u = 1.5$
- Mass balance $\sim10^{-15}$
- Pressure drop vs $\Delta p = 12\mu\bar u L/h^2$

The pressure-drop check is the sharpest of these, because it exercises the
pressure solve, the boundary conditions, and the viscous term together. If
$\Delta p$ is right to a fraction of a percent, your channel solver works.

**Do not skip this step when you build your own geometry.** Every new geometry
should start as a straight channel that you verify, and only then get its
obstacle/bend/contraction added. Debugging a serpentine mixer from scratch is
miserable; debugging one change to a working channel is easy.

---

## 8.3 Solid obstacles: the masking approach

The simplest way to put a solid body in a Cartesian grid is **masking** (a crude
immersed-boundary method):

1. Define a logical array `solid(i,j)` marking cells inside the body.
2. Zero the velocity on any face touching a solid cell.
3. **Put the mask into the pressure operator too.**

```matlab
% zero every u-face that borders a solid cell
uMask = false(nx+1, ny);
uMask(1:end-1,:) = uMask(1:end-1,:) | solid;
uMask(2:end,:)   = uMask(2:end,:)   | solid;
u(uMask) = 0;
```

### Step 3 is not optional, and step 3 is the one people skip

The natural assumption is that zeroing velocities is enough — "the zero velocity
constrains the pressure correctly." It does not. **Measured in this code:
masking velocities alone gave `max|div u| = 416` and a 9.7% mass-balance error;
putting the mask into the Poisson operator brought both to $10^{-14}$.**

The reason is simple once stated. The projection makes the flow divergence-free
by construction — but only with respect to the operator it solved. If that
operator thinks fluid can flow through the body, it produces a $\phi$ whose
gradient drives flow into the solid; you then zero those faces afterwards, and
the cells next to the body are left with whatever divergence that zeroing
created. **You cannot make a flow divergence-free around a body the pressure
solve does not know about.**

The fix is that a solid face and a Neumann wall are the *same thing* — both are
zero-flux — so the operator just drops that face:

```matlab
elseif mask(i+1,j)
    % solid neighbour: no flux, contribute nothing (exactly like a Neumann wall)
else
    add(k, k+1, cx);  d = d - cx;
end
```

See `poisson2d.m`, which takes an optional mask argument for this purpose.

**Honest limitations, because this method is still easy to over-trust:**

- The boundary is **staircased** — a circle becomes a pixelated blob. Forces on
  the body are consequently first-order accurate at best, and the drag
  coefficient you compute will be a few percent off even on a fine grid.
- No-slip is enforced at the *face*, not at the true surface, so there is a
  half-cell geometric error.
- You need enough cells across the body: **at least 20** across the smallest
  feature, more if you care about shear stress on it.

For microfluidics this is usually acceptable — most microchannel geometries are
rectilinear anyway (they're made by photolithography, so the walls really *are*
axis-aligned). For a cylinder in cross-flow where you need accurate drag, use a
proper cut-cell or immersed-boundary method, or a body-fitted mesh.

`code/lib/make_mask.m` builds masks from simple primitives (rectangles, circles,
and function handles) so you can compose geometry without writing index logic.

---

## 8.4 Three worked geometries

### Contraction / expansion (`channel.m`, `'contraction'`)

A channel that narrows from $h_1$ to $h_2$. Conservation of mass gives
$\bar u_2 = \bar u_1 h_1/h_2$, so this is a self-checking case.

**Microfluidic relevance:** contractions are where extensional flow lives.
Polymers and cells stretch here; it's also where clogging starts.

At $Re \ll 1$ there is **no separation** at the expansion — the flow follows the
geometry exactly. At $Re > 10$ you'd get a recirculation bubble. Running both is
a good way to see what $Re$ actually buys.

### Cylinder in a channel (`obstacle.m`)

The classic. At $Re \ll 1$ the streamlines are **fore-aft symmetric** — the
picture upstream is a mirror image of downstream. That symmetry is the visual
signature of Stokes reversibility (Tutorial 1) and it is a *strong* test: any
visible asymmetry at $Re = 0.01$ means your convective term has a bug or is
being fed a wrong velocity.

At $Re = 40$ you get standing recirculation bubbles behind the cylinder. The
transition happens around $Re \approx 5$.

### Serpentine (`serpentine.m`)

A channel that snakes back and forth. The workhorse geometry of passive
micromixers, and the input to Tutorial 9.

**What to look for:** at each bend, the streamlines on the inside of the turn
are compressed and move faster. At high $Re$ this generates **Dean vortices**
(secondary flow in the cross-section) which mix well — but Dean flow is a 3D
effect that a 2D simulation *cannot* capture. At $Re \ll 1$ there is no Dean
flow anyway, so the 2D limitation is not binding, but you should know it's there.

The Dean number is $De = Re\sqrt{h/(2R)}$; below $De \approx 10$ there is
essentially no secondary flow.

---

## 8.5 Choosing a resolution

Three constraints, and you must satisfy all of them:

**1. Resolve the geometry.** ≥ 20 cells across the smallest feature.

**2. Resolve the velocity gradients.** ≥ 20 cells across the channel to get the
parabolic profile and hence wall shear stress right. (10 gets the flow rate but
not the shear.)

**3. Control numerical diffusion** — usually the binding one if you will
transport a scalar. From Tutorial 5, $Pe_{cell} = u\Delta x/D$, and you want
$D_{num}/D < 0.1$.

Constraint 3 is brutal. For the standard microchannel it demands
$\Delta x < 0.1$ µm with upwind — 1000 cells across a channel. With van Leer
(Tutorial 5 measured ~750× less numerical diffusion) you can get away with
50–100 cells. **This is the reason the mixing tutorial uses a TVD scheme; on a
grid you can actually afford, nothing else is honest.**

Budget: a 2D run of 400×100 cells is ~40k unknowns — seconds per hundred steps
in MATLAB. 2000×200 is 400k unknowns, and you will be waiting. Start coarse,
establish that the physics looks right, then refine and check that your answer
doesn't move.

---

## 8.6 A workflow that works

1. **Sketch the geometry**, mark inlet/outlet/wall on each boundary.
2. **Non-dimensionalize** and compute $Re$, $Pe$, $Pe_{cell}$ (Tutorial 1's
   `scaling`). Decide the grid from §8.5.
3. **Run a straight channel** of the same width and grid. Verify Poiseuille,
   mass balance, pressure drop. *Do not skip.*
4. **Add the geometry.** Re-check mass balance and `divMax`.
5. **Grid-convergence check.** Refine 2×. If your answer moves more than a few
   percent, it wasn't converged. Report the number you get *and* the number you
   got on the coarse grid.
6. **Sanity-check against theory.** Pressure drop vs the analytic estimate,
   symmetry at low $Re$, flow rate splits at junctions.
7. **Only then** add scalar transport (Tutorial 9) or two-phase (Tutorial 11).

Step 5 is the one people skip and the one that invalidates results. A CFD number
without a grid-convergence check is a guess with error bars you didn't measure.

---

## Exercises

**8.1** Run `channel`. Confirm the outlet profile is parabolic, mass balance is
$\sim10^{-15}$, and pressure drop matches $12\mu\bar u L/h^2$ to under 1%.

**8.2** Break the outlet consistency deliberately: replace the Dirichlet-face
correction in `ns_solver` with `u(end,:) = u(end-1,:)` and rerun `channel`.
Watch `divMax` jump from ~1e-14 to ~1. Then localize it — plot `div` as an image
and confirm the error lives entirely in the last column of cells. *(Do this once
so you recognize the symptom: a conservation error confined to one boundary is
always a boundary-consistency bug, never a stability or resolution problem.)*

**8.2b** Same experiment for the mask: pass `[]` instead of `opt.solid` to
`poisson2d` inside `ns_solver`, keep the velocity masking, and rerun the
contraction. `divMax` should go to ~1e2. Plot where it lives.

**8.3** Run `obstacle` at $Re = 0.01$. Measure the fore-aft asymmetry
quantitatively: compare $u(x_c - d)$ with $u(x_c + d)$ for several $d$. How
symmetric is it? Now run $Re = 40$ and repeat.

**8.4** Run `channel` with `'contraction'` from 100 µm to 50 µm. Verify
$\bar u$ doubles. Compute the extra pressure drop above the sum of the two
straight sections — that's the *minor loss* of the contraction.

**8.5** Compute the drag on the cylinder in `obstacle` by summing pressure and
viscous forces over the masked faces. Compare to Stokes drag. Then refine the
grid 2× and see how much the answer moves. *(This is where the staircasing
limitation of §8.3 shows up. Report both numbers.)*

**8.6 (your own device)** Build the geometry from your Exercise 1.5 device.
Follow the §8.6 workflow all the way through step 6. Save the resulting `S`
struct — Tutorial 9 will transport a scalar on it.

---

**Next:** [Tutorial 9 — Mixing and dispersion](09-mixing-and-dispersion.md).
You have a flow field. Now find out how badly it mixes, which is the central
question of microfluidics.
