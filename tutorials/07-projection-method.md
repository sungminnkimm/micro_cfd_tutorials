# Tutorial 7 — The Projection Method

**Goal:** Write a complete incompressible Navier–Stokes solver. This is the
centrepiece of the course; everything in Part I was preparation and everything
in Part III is application.

**Code:** `code/lib/ns_solver.m`, `code/ch07/cavity.m`, `code/ch07/cavity_validate.m`

---

## 7.1 The problem with pressure

$$\frac{\partial\mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu\nabla^2\mathbf{u} + \mathbf{f}$$
$$\nabla\cdot\mathbf{u} = 0$$

Three unknowns ($u, v, p$), three equations. But as Tutorial 6 noted, there is no
evolution equation for $p$ — the second equation is a *constraint*, not a
transport equation. You cannot march pressure forward in time because there is
nothing to march.

**Chorin's insight (1968):** don't try. Split the step in two.

1. **Predictor.** Ignore pressure. Advance velocity using only convection,
   diffusion, and body forces. Call the result $\mathbf{u}^*$. It will not be
   divergence-free.
2. **Projector.** Find the pressure field that, applied as a correction,
   removes exactly the divergence you just created.

---

## 7.2 The mathematics: Helmholtz–Hodge

The reason this works is a theorem. **Any** vector field can be decomposed
uniquely into a divergence-free part and a gradient:

$$\mathbf{u}^* = \underbrace{\mathbf{u}^{n+1}}_{\nabla\cdot\ =\ 0} + \nabla\phi$$

So the divergence-free part is what remains after you subtract off a gradient —
$\mathbf{u}^{n+1}$ is the **projection** of $\mathbf{u}^*$ onto the space of
divergence-free fields. Hence the name.

To find $\phi$, take the divergence of both sides and use
$\nabla\cdot\mathbf{u}^{n+1} = 0$:

$$\nabla\cdot\mathbf{u}^* = \nabla^2\phi$$

**A Poisson equation.** The same one from Tutorial 3, which you have already
built, verified by MMS, and factorized. That is the whole trick: the hardest
part of incompressible flow reduces to a problem you already solved.

### The algorithm

$$\text{(1)}\quad \mathbf{u}^* = \mathbf{u}^n + \Delta t\left[-(\mathbf{u}\cdot\nabla)\mathbf{u} + \nu\nabla^2\mathbf{u} + \mathbf{f}\right]^n$$
$$\text{(2)}\quad \nabla^2\phi = \frac{\nabla\cdot\mathbf{u}^*}{\Delta t}$$
$$\text{(3)}\quad \mathbf{u}^{n+1} = \mathbf{u}^* - \Delta t\,\nabla\phi$$

with $p = \rho\phi$ (to first order; see §7.6). Three steps, one Poisson solve.

**The check that proves it worked:** after step 3, compute
$\nabla\cdot\mathbf{u}^{n+1}$. On a staggered grid it should be $\sim10^{-14}$.
Not "small" — *machine precision*. If it isn't, something is wrong, and the
solver prints this number every so often for exactly that reason.

---

## 7.3 Boundary conditions for $\phi$ — the subtle part

What boundary condition does the Poisson equation for $\phi$ take? This trips up
almost everyone, and getting it wrong produces a solver that *almost* works.

Take the normal component of step (3) at a wall:

$$u^{n+1}_n = u^*_n - \Delta t\frac{\partial\phi}{\partial n}$$

At a solid wall, $u^{n+1}_n = 0$ (no penetration). If you also enforce
$u^*_n = 0$ when building the predictor — which you should, and which the code
below does — then

$$\boxed{\frac{\partial\phi}{\partial n} = 0 \quad\text{at every solid wall}}$$

**Homogeneous Neumann.** Which means, for a fully enclosed domain like the
lid-driven cavity: **all-Neumann, hence singular**, exactly the case Tutorial 3
warned about. You pin one cell. `poisson2d` already does this automatically when
it sees four `'N'` walls.

And the compatibility condition $\int f\,dV = \oint g\,dS$ becomes
$\int \nabla\cdot\mathbf{u}^*\,dV = 0$, i.e. **net flux through the boundary must
be zero**. For a closed box that's automatic. For the channel of Tutorial 8, it
means inflow must exactly equal outflow — and if it doesn't, the pressure solve
fails. That's the practical bite of the theory from Tutorial 3.

At an **outflow** boundary the story differs: you specify pressure instead, so
$\phi$ gets a Dirichlet condition there, and the singularity goes away. Tutorial 8.

---

## 7.4 The convective term on a staggered grid

This is the fiddliest code in the course. We use the **divergence (conservative)
form**, which conserves momentum exactly:

$$(\mathbf{u}\cdot\nabla)u = \frac{\partial(u^2)}{\partial x} + \frac{\partial(uv)}{\partial y}$$

To evaluate this at a $u$-face you need $u^2$ at **cell centres** and $uv$ at
**cell corners**:

```matlab
% u^2 at cell centres: average the two u-faces bounding each cell
uc = 0.5*(ue(1:end-1,:) + ue(2:end,:));        % nx x ny

% uv at corners: average u vertically, v horizontally
u_cor = 0.5*(ue(:,1:end-1) + ue(:,2:end));     % (nx+1) x (ny+1)
v_cor = 0.5*(ve(1:end-1,:) + ve(2:end,:));     % (nx+1) x (ny+1)

% now differentiate back onto the u-faces
conv_u = (uc(2:end,:).^2 - uc(1:end-1,:).^2)/dx ...
       + (u_cor(2:end-1,2:end).*v_cor(2:end-1,2:end) ...
        - u_cor(2:end-1,1:end-1).*v_cor(2:end-1,1:end-1))/dy;
```

where `ue`, `ve` are the ghost-padded arrays. **Draw the grid and check every
index.** The rule that makes it tractable: *interpolate to where the product
lives, multiply there, then difference back to where you need the result.* Never
interpolate a product.

At the low $Re$ of microfluidics this term is nearly negligible — which is
convenient, because it means a bug here is easy to miss. Test your solver at
$Re = 100$ (lid-driven cavity, §7.7) where the term *does* matter, even if you
only ever run it at $Re = 0.1$.

---

## 7.5 Timestep

Three limits, all of which you have now derived:

```matlab
dt_conv = dx / max(abs(u(:)));                     % CFL          (Tutorial 5)
dt_visc = 0.5 / nu / (1/dx^2 + 1/dy^2);            % viscous      (Tutorial 4)
dt = safety * min(dt_conv, dt_visc);               % safety ~ 0.5
```

**At the microscale, `dt_visc` almost always wins.** With $\nu = 10^{-6}$ m²/s
and $\Delta x = 1$ µm, $\Delta t \le 2.5\times10^{-7}$ s while the CFL limit at
1 mm/s is $10^{-3}$ s — a factor of 4000. This is the practical face of $Re\ll1$
and it is why:

- Explicit viscous treatment is expensive at the microscale.
- Tutorial 9 solves the flow **once** to steady state and then freezes it.

If you need long simulations, treat the viscous term implicitly (Crank–Nicolson,
Tutorial 4) and only the convection explicitly. The solver below supports this
via `opt.implicitViscous`.

---

## 7.6 Accuracy caveats worth knowing

**$\phi$ is not exactly $p/\rho$.** Chorin's scheme is first-order in time, and
the relationship is $p = \rho\phi + O(\Delta t)$ — with a viscous correction
$p = \rho\phi - \tfrac{\mu\Delta t}{2}\nabla^2\phi$ in the second-order variants.
For steady flows this doesn't matter (everything converges to the right answer).
For unsteady pressure histories, it does.

**The numerical boundary layer.** The $\partial\phi/\partial n = 0$ condition is
not exactly what the continuous problem wants; it induces a thin $O(\sqrt{\nu\Delta t})$
error layer at walls. It's small and it doesn't pollute the interior, but it does
mean the *pressure* converges at a lower rate than the velocity near walls.

**Second-order variants** (van Kan / incremental pressure-correction) fix both by
carrying the previous pressure:

$$\mathbf{u}^* = \mathbf{u}^n + \Delta t\left[\dots - \tfrac{1}{\rho}\nabla p^n\right], \qquad p^{n+1} = p^n + \rho\phi$$

One extra line, second-order in time. `ns_solver` uses this by default
(`opt.incremental = true`).

---

## 7.7 Validation: the lid-driven cavity

The standard benchmark for incompressible solvers. A square box, three
stationary walls, a lid sliding at constant $U$. The published reference is
**Ghia, Ghia & Shin (1982)**, tabulated centreline velocities at several
Reynolds numbers, and every new solver gets checked against it.

`code/ch07/cavity_validate.m` runs $Re = 100$ and $Re = 400$ and compares
against the tabulated Ghia data (embedded in the file, so no download needed).

**Actual measured result** (this is what the code in this repo produces, not an
aspiration):

```
    Re    N     err(u)    err(v)     max|div u|
   ---------------------------------------------
   100    64    0.0038    0.0085      7.11e-15
   100   128    0.0049    0.0091      1.37e-14
   400    64    0.0073    0.1435      7.11e-15
   400   128    0.0020    0.1494      1.33e-14
```

Errors are absolute, with the lid speed equal to 1 — so $Re=100$ agrees with the
published data to better than 1% of the lid speed, and `max|div u|` is at
round-off. That is a pass.

### The Re = 400 anomaly, and why it is worth showing you

`err(v)` at $Re=400$ is 0.14, and **it does not improve when the grid is
refined** (0.1435 → 0.1494), while `err(u)` on the same runs converges cleanly
(0.0073 → 0.0020).

That asymmetry is diagnostic, and reasoning it through is a better lesson than a
clean table would be. A genuine under-resolution error *shrinks* when you refine.
An error that stays flat under refinement is not a discretization error at all.
And since $u$ and $v$ come from the *same* velocity field, a solver bug that
wrecked $v$ while leaving $u$ converging at second order is hard to construct.

The likely culprit is the reference data — the $Re=400$ $v$ row embedded in
`ghia_data.m` was transcribed from memory and is **not** verified against the
original paper. It is flagged as such in the file.

**Check it against Table II of Ghia et al. before using it to judge your own
code.** The general point matters more than this particular row: benchmark data
deserves the same scepticism as your solver. When code and reference disagree,
"the code is wrong" is a hypothesis, not a conclusion — and the
does-it-converge-under-refinement test is how you tell the two apart.

**Why validate at $Re = 100$ when microfluidics is $Re = 0.1$?** Because at
$Re = 0.1$ the convective term is invisible — you could delete it entirely and
the cavity result would barely change. A test that cannot fail proves nothing.
Validating where the physics is *hard* is what tells you the code is right when
you later run it where the physics is easy.

### Diagnostics to watch

The solver prints, every `opt.report` steps:

| Quantity | Healthy value | If wrong |
|---|---|---|
| `max\|div u\|` | $10^{-14}$ | Poisson solve or BC bug |
| `max\|u\|` | O(lid speed) | growing → unstable |
| `dt` | steady | shrinking → something is accelerating |
| `d(KE)/dt` | → 0 at steady state | not converging |

**`max|div u|` is the single most valuable number.** If it is at round-off, your
projection is working, and any remaining error is in the predictor. That splits
the debugging problem in half.

---

## 7.8 The solver

`code/lib/ns_solver.m` is ~180 lines. Its structure:

```
setup:   grid, factorize the pressure Poisson operator ONCE
loop:
   1. apply velocity BCs (build ghost-padded arrays)
   2. compute convection + diffusion -> u*, v*
   3. apply BCs to u*, v*  (normal components must match u^{n+1})
   4. rhs = div(u*)/dt
   5. phi = poisson_solve(...)        <- the precomputed factorization
   6. u = u* - dt*grad(phi)
   7. check div(u), report
```

Read the file alongside this list. The one thing to note is **step 3**: the
normal component of $\mathbf{u}^*$ at a wall must be set to the same value
$\mathbf{u}^{n+1}$ will have. If you skip it, the $\partial\phi/\partial n = 0$
derivation in §7.3 is invalid, and you get a solver that produces plausible
flows with a slow mass leak. That's the bug that takes a week to find; set the
normal components and it never happens.

### Using it

```matlab
opt = ns_defaults();
opt.nx = 64;  opt.ny = 64;
opt.Lx = 1e-3;  opt.Ly = 1e-3;
opt.nu = 1e-6;  opt.rho = 998;
opt.bc.top = struct('type','moving','u',1e-3);   % lid
opt.tEnd = 2.0;

S = ns_solver(opt);
quiver(S.G.xp, S.G.yp, S.uc.', S.vc.');    % note the transposes
```

---

## Exercises

**7.1** Run `cavity` at $Re = 100$. Confirm `max|div u|` is at round-off.
Deliberately comment out the projection step (6) and watch it grow.

**7.2** Run `cavity_validate`. Compare with Ghia at $Re = 100$ on 64² and 128².
Does refining help? Quantify.

**7.3 (the one that finds bugs)** Set $Re = 0$ by making `nu` enormous, so the
convective term is negligible. Then set the lid velocity to zero and give the
flow a random initial condition. It must decay to zero. Anything that persists
is a spurious mode.

**7.4** Time one step, broken down into predictor / Poisson / corrector. Which
dominates? Now re-run with the factorization moved *inside* the loop and measure
the slowdown. *(This is the §3.4 claim, verified on the real solver.)*

**7.5** Reduce to a microfluidic case: $L = 100$ µm, lid at 1 mm/s, water. Compute
$Re$. Compare the streamline pattern to $Re = 100$. Where has the vortex centre
moved, and why? *(The answer — toward the geometric centre — is the visual
signature of Stokes flow's reversibility.)*

**7.6** Switch on `opt.implicitViscous` and compare the allowed timestep and total
runtime against the explicit version for the microfluidic case in 7.5.

**7.7** Verify the solver's *temporal* order using MMS: pick an exact unsteady
solution (Taylor–Green vortex,
$u = -\cos x\sin y\,e^{-2\nu t}$, $v = \sin x\cos y\,e^{-2\nu t}$, which solves NS
exactly on a periodic domain), and confirm second order with
`opt.incremental = true` and first order without. *(Tutorial 13 gives the full
recipe if you want to do this properly.)*

---

**Next:** [Tutorial 8 — Building your own geometry](08-microchannel-geometry.md).
You have a solver. Now learn to point it at a channel you designed, with inlets,
outlets, and obstacles.
