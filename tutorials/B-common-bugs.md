# Appendix B — Common Bugs and Their Symptoms

A lookup table from *what you see* to *what's wrong*. Ordered by how often each
one actually happens.

---

## B.1 The master diagnostic table

| Symptom | Cause | Where |
|---|---|---|
| `NaN` appears in **one step** | Division by zero: unguarded limiter ratio, or $1/\|\nabla\phi\|$ | T5.5, T11 |
| Blows up over **~10 steps**, cell-to-cell sawtooth | Timestep exceeds a stability limit | T4.2 |
| Cell-to-cell **wiggles that don't grow** | $Pe_{cell}>2$ with central differencing | T5.4 |
| Smooth, plausible, **wrong** | Boundary condition. Almost always | T2.3, T3.3 |
| Converges at **order 1** instead of 2 | Boundary condition, or a first-order stencil leaking in | T2.3 |
| Order **drops on the finest grids** | Round-off floor — expected, stop refining | T2.3 |
| `divMax` **not at round-off** | Pressure operator inconsistent with the corrector | T8.1 |
| Mass **not conserved** | Outflow BC, or level-set reinit drift | T8.1, T11.3 |
| Result **changes with grid** | Not converged. A measurement, not a bug | T13.4 |
| Diverges at a **safe timestep**, `divMax` fine | $Re_{cell}=u\Delta x/\nu>2$ — use `opt.stokes` | §B.4 |
| Mixing looks **too good** | Numerical diffusion | T5.6 |
| Interface **wobbles** | Spurious currents | T11.4 |
| Droplet **dissolves** | Level-set mass loss; enable the area correction | T11.3 |
| Everything is **transposed** | `meshgrid` instead of `ndgrid`, or swapped `kron` args | T3.5 |
| Works on square grids, **fails on rectangular** | Same as above — and that's how you find it | T3.2 |

---

## B.2 Distinguishing the three failure classes

This is the single most useful diagnostic skill, because the three look
superficially similar and have completely different fixes.

**Instability** — cell-to-cell sawtooth that **grows geometrically**, then
`Inf`.
→ Reduce the timestep. Check *every* limit (CFL, viscous, capillary), not just
the one you remember.

**Wiggles** — cell-to-cell oscillation that is **bounded and steady**, near
sharp gradients only.
→ Not instability. $Pe_{cell}>2$. Refine, or switch to upwind/TVD.

**Boundary bug** — solution is **smooth and plausible** everywhere, correct in
the interior, wrong near one edge; convergence rate is 1 instead of 2.
→ The hardest to spot by eye and the easiest to catch by an order-of-accuracy
test. This is the argument for MMS.

```matlab
% Quick classifier
saw = max(abs(diff(c,2,1)),[],'all');   % cell-to-cell curvature
if ~isfinite(saw),        disp('instability');
elseif saw > 0.5*range(c(:)), disp('wiggles or instability - check growth');
else,                     disp('smooth: suspect boundaries');
end
```

---

## B.3 Staggered-grid specific

These account for most of the pain in Tutorials 6–12.

| Bug | Symptom | Check |
|---|---|---|
| Wrong array size (`nx` vs `nx+1`) | Immediate size error, or silent broadcasting | `size(u)` vs `G.size_u` |
| Ghost cells on the wrong boundaries | Wrong profile at two walls only | Tutorial 6.3's table |
| `u_ghost = -u_in` instead of `2*u_w - u_in` | No-slip works, moving wall doesn't | Test with a moving lid |
| Interpolating a product instead of multiplying interpolants | Subtly wrong convection; only visible at higher $Re$ | Validate at $Re=100$, not $Re=0.1$ |
| Half-cell offset between $u$ and $p$ | Order drops to 1 | MMS |

**The rule that prevents most of these:** *interpolate to where the product
lives, multiply there, then difference back.* Never interpolate a product.

**And: validate where the physics is hard.** At $Re=0.1$ you could delete the
convective term entirely and the answer would barely change — so a test at
$Re=0.1$ cannot detect a bug in it. This is why Tutorial 7 validates against
Ghia at $Re=100$ and $400$ even though the course runs at $Re\ll1$.

---

## B.4 Three bugs from building this course

All three shipped plausible-looking output. None crashed. Two gave a *correct*
answer to the first quantity anyone would check.

### 1. Solid mask missing from the pressure operator

**Symptom:** contraction geometry ran fine, throat velocity exactly 2.001×
upstream as mass conservation demands — but `divMax` was **416**.

**Cause:** the mask was applied to velocities after the projection, but not to
the Poisson operator. The projection was making the flow divergence-free with
respect to a domain containing no obstacle.

**Fix:** a solid face *is* a zero-flux face, so the operator simply drops it —
identical treatment to a Neumann wall. `divMax` went to $7.8\times10^{-14}$.

**Lesson:** you cannot make a flow divergence-free around a body the pressure
solve doesn't know about.

### 2. Outflow face not pressure-corrected

**Symptom:** after fixing (1), `divMax` was **1.5**. Better, not round-off.

**Cause:** the outlet velocity was re-extrapolated *after* the projection
instead of being corrected with the Dirichlet-$\phi$ ghost the operator had
assumed. The error lived entirely in the last column of cells.

**Fix:** correct the boundary face with the matching gradient,
`u(end,:) = us(end,:) - dt*(0 - phi(end,:))*2/dx`. → $4.2\times10^{-14}$.

**Lesson:** "1.5 is not $10^{-14}$" must be treated as meaningful. A
conservation error confined to one boundary is *always* a boundary-consistency
bug, never a resolution problem.

### 3. Curvature sign

**Symptom:** static droplet gave a pressure jump of **−169 Pa** where theory
says **+125 Pa**.

**Cause:** with $\phi>0$ inside, $\nabla\phi$ points *inward*, so the standard
formula returns $\kappa = -1/R$. Surface tension was pushing the droplet apart.

**Fix:** negate, so $\kappa = +1/R$ for a droplet and the CSF force points
inward. Spurious currents fell by ~30 orders of magnitude and Δp came to 1.0% of
$\sigma/R$.

**Lesson:** there was an exact answer available ($\sigma/R$), and checking it
rather than eyeballing the picture is what caught it.

### 4. Pressure drift at an outlet — and a wrong diagnosis on the way

Worth reading in full, because the *misdiagnosis* is as instructive as the bug.

**Symptom:** the serpentine channel diverged at step 167, with a timestep at
*half* the viscous stability limit. `max|div u|` stayed at $10^{-13}$ the entire
time — the projection was working flawlessly right up to the blow-up — while
`max|u|` grew ×5 every 25 steps.

**First (wrong) hypothesis:** cell-Reynolds instability. Central differencing of
convection is unconditionally unstable alone (Tutorial 5.2) and is stabilized
only by viscosity, requiring $Re_{cell} = |u|\Delta x/\nu \lesssim 2$. A quick
mental estimate gave $Re_{cell}\approx7.5$ — comfortably past the limit, and a
satisfying story.

Setting `opt.stokes = true` (dropping convection entirely, which $Re\ll1$
justifies anyway) made it survive **36× longer** — 6005 steps instead of 167.
Encouraging. And still wrong: it then diverged anyway.

The estimate was off by a factor of 1000:

$$Re_{cell} = \frac{1.5\times10^{-3}\times5\times10^{-6}}{1.0\times10^{-6}} = 0.0075$$

nowhere near 2. The improvement from `stokes` was real but incidental — removing
a term removed some of the energy feeding the true instability, which delayed it
without curing it.

**Actual cause:** printing *where* the blow-up lived settled it in one line — at
`u(401,49)`, the **outlet face**, with `max|p| = 5.9e154`.

In the incremental (van Kan) scheme the pressure accumulates,
$p^{n+1} = p^n + \rho\phi$. Setting $\phi=0$ at a pressure outlet pins the
**correction** to zero — not the **pressure**. So whatever $p$ has drifted to at
the outlet simply stays there and keeps growing. The straight channel hid it by
reaching steady state after 1936 steps, before the drift became visible.

**Fix:** to make $p^{n+1}=0$ at the outlet, give $\phi$ the wall value
$\phi_{wall} = -p_{wall}/\rho$ rather than 0 — and use the *same* value in the
face correction, or operator and corrector disagree (bug #2 again, one level up).

**Two lessons.**

1. **A healthy `divMax` does not mean a healthy solver.** It certifies the
   projection and nothing else. Every other term can be diverging beneath it.
2. **Locate before you theorize.** `[~,k]=max(abs(u(:))); ind2sub(...)` took
   seconds and pointed straight at the outlet. An hour went into a plausible
   theory that a single `fprintf` of the actual number would have killed
   immediately. When you have a hypothesis with a number in it, *print the
   number*.

### 5. Inflow profile orientation

**Symptom:** `Number of elements must not change` from `reshape`, in the
two-phase solver only.

**Cause:** an inflow profile is naturally built as a **column** (from `G.yu`),
but the slice it fills, `us(1,:)`, is a **row**. Writing
`profile .* ones(size(inner))` then broadcasts column against row and silently
produces an $n_y\times n_y$ **matrix**.

**Lesson:** MATLAB's implicit expansion turns a shape mismatch into a silently
larger array instead of an error. When a reshape complains about element counts,
suspect broadcasting upstream — the bug is rarely at the reshape itself.

### 6. Level-set reinitialization dissolving the droplet

**Symptom:** area error exactly $-1.0$ — the droplet vanished completely.

**Cause:** each reinitialization drifts the zero contour by ~0.1% of the area.
Over 3,000 calls, $0.999^{3000}\approx0$.

**Fix:** global area correction after each reinit. Error went from $-1.0$ to
$-1.5\times10^{-6}$.

**Lesson:** a per-call error of 0.1% is not small when you make it 3,000 times.
Always check the *cumulative* budget of a repeated operation.

---

## B.5 Fast sanity checks

Run these on any new solver, in this order. Each takes seconds.

```matlab
% 1. Zero everything -> nothing should happen
% 2. Uniform flow -> stays uniform
% 3. Symmetric geometry + symmetric BCs -> symmetric solution
% 4. Reverse all BCs -> reversed solution (Stokes only)
% 5. Refine 2x -> error drops 4x
% 6. div(u) after projection ~ 1e-14
% 7. sum(c) with closed walls -> constant
```

Checks 1–4 cost nothing and catch gross errors. Checks 5–7 are the ones that
catch the subtle bugs, and they are the ones people skip.

---

## B.6 When you are truly stuck

1. **Shrink it.** 8×8 grid. Print the whole matrix. A bug in a 400×100 droplet
   run is unfindable; the same bug in an 8×8 static case is obvious.
2. **Remove terms.** Turn off convection. Turn off surface tension. Find the
   smallest configuration that still fails.
3. **Check the operator, not the solution.** `full(A)` on a tiny grid, `spy(A)`,
   check row sums (a pure-Neumann Laplacian row must sum to zero).
4. **Compare against something exact.** MMS lets you manufacture an exact answer
   for *any* PDE — you are never without a reference.
5. **Question the reference.** Tutorial 7's Ghia comparison found a discrepancy
   that did not shrink under refinement, which is the signature of bad
   *reference data*, not a bad solution. When code and reference disagree, "the
   code is wrong" is a hypothesis, not a conclusion.
