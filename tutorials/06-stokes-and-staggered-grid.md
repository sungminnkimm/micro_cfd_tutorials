# Tutorial 6 — Stokes Flow and the Staggered Grid

**Goal:** Understand why storing pressure and velocity at the same place breaks,
adopt the staggered (MAC) grid that fixes it, and solve your first real flow
problem with an analytic answer to check against.

**Code:** `code/ch06/checkerboard_demo.m`, `code/ch06/mac_grid.m`, `code/ch06/poiseuille.m`

---

## 6.1 The equations

Incompressible flow, in the low-$Re$ limit that Tutorial 1 justified:

$$\nabla p = \mu\nabla^2\mathbf{u}, \qquad \nabla\cdot\mathbf{u} = 0$$

Notice something odd about that pair. There is **no time derivative of pressure**
and **no evolution equation for pressure at all**. Pressure is not a
thermodynamic state variable here — it is a **Lagrange multiplier** that
enforces $\nabla\cdot\mathbf{u} = 0$. It takes whatever value is needed,
instantaneously, everywhere, to keep the flow divergence-free.

This is the central structural difficulty of incompressible CFD, and everything
in Tutorials 6 and 7 is machinery for coping with it.

---

## 6.2 The checkerboard problem

Suppose you store $u$, $v$, and $p$ all at cell centres — a **colocated** grid,
the obvious choice — and use central differences for the pressure gradient:

$$\frac{\partial p}{\partial x}\bigg|_i = \frac{p_{i+1} - p_{i-1}}{2\Delta x}$$

Now consider this pressure field, alternating cell to cell:

```
   p:   1   0   1   0   1   0   1   0
```

Its central-difference gradient at every node is $(0 - 0)/2\Delta x$ or
$(1-1)/2\Delta x$ — **exactly zero everywhere**.

The momentum equation cannot see this field. It is invisible: a nonzero pressure
distribution that exerts no force. It sits in the null space of your discrete
gradient operator, so the solver can add any amount of it to the solution and
nothing objects. In practice it grows from round-off until your pressure field is
a wild checkerboard while the velocity field looks deceptively fine.

The velocity divergence has the mirror problem: the same central stencil applied
to a checkerboard *velocity* field reports zero divergence, so a physically
impossible flow passes the incompressibility check.

`code/ch06/checkerboard_demo.m` constructs the mode and confirms the discrete
gradient of it is zero to machine precision — worth running once so you believe
the problem is real and not folklore.

### Two ways out

1. **Rhie–Chow interpolation** — keep the colocated grid, add a carefully
   designed pressure-smoothing term to the face velocities. This is what most
   commercial codes (Fluent, OpenFOAM) do, because unstructured meshes make
   staggering awkward. It works, but it introduces a tunable damping and it is
   fiddly to derive.
2. **Staggered grid** — store the variables at different locations so the
   problem never arises. Simpler, exactly conservative, and the natural choice
   on a structured Cartesian grid.

**This course uses staggering.** Harlow & Welch introduced it in 1965 for the
Marker-and-Cell method, hence "MAC grid."

---

## 6.3 The MAC grid

The idea: put each variable where its derivative is naturally needed.

- **Pressure** at cell **centres** ($\times$)
- **$u$** (x-velocity) at the **vertical faces** ($\rightarrow$)
- **$v$** (y-velocity) at the **horizontal faces** ($\uparrow$)

```
        v(i,j+1)
           ↑
    ┌──────┴──────┐
    │             │
 →  │      ×      │  →     u(i,j) ... p(i,j) ... u(i+1,j)
u(i,j)   p(i,j)    u(i+1,j)
    │             │
    └──────┬──────┘
           ↑
        v(i,j)
```

Now look what happens to the two problematic operators.

**Pressure gradient**, needed at the $u$-face between cells $i-1$ and $i$:

$$\frac{\partial p}{\partial x}\bigg|_{u(i,j)} = \frac{p_{i,j} - p_{i-1,j}}{\Delta x}$$

A **compact two-point difference** across one cell — not a wide four-cell-apart
one. It is second-order accurate *at the face where it is evaluated*, and a
checkerboard pressure now produces a large alternating gradient rather than
zero. The null space is gone.

**Divergence**, needed at the cell centre:

$$(\nabla\cdot\mathbf{u})_{i,j} = \frac{u_{i+1,j} - u_{i,j}}{\Delta x} + \frac{v_{i,j+1} - v_{i,j}}{\Delta y}$$

Also compact, and — this is the part that matters physically — it is **exactly**
the net volume flux through the four faces of cell $(i,j)$, divided by the cell
volume. Mass conservation is not approximated. It is enforced to machine
precision, cell by cell.

That exactness is why staggered grids are worth the indexing pain. In Tutorial 7
you will check `max(abs(div))` after every projection and expect $10^{-14}$; on a
colocated grid you would have to settle for "small."

### Array sizes — write these on a sticky note

For a domain of $n_x \times n_y$ **cells**:

| Variable | Size | Located at |
|---|---|---|
| `p` | `nx × ny` | cell centres |
| `u` | `(nx+1) × ny` | vertical faces, **including both x-boundaries** |
| `v` | `nx × (ny+1)` | horizontal faces, **including both y-boundaries** |

The `+1` is the entire source of staggered-grid bugs. `u` has one more column
than `p` because a row of $n_x$ cells has $n_x+1$ faces — think of fenceposts.

Coordinates:

```matlab
xp = (dx/2 : dx : Lx-dx/2);     % nx    pressure/cell centres
yp = (dy/2 : dy : Ly-dy/2);     % ny
xu = (0    : dx : Lx);          % nx+1  u-faces
yu = yp;                        % ny    u sits at cell-centre height
xv = xp;                        % nx
yv = (0    : dy : Ly);          % ny+1  v-faces
```

**A boundary condition on $u$ at $x=0$ is now a Dirichlet condition on an actual
unknown**, `u(1,:)` — no ghost cell, no interpolation. That is a real
simplification. But $u$ at the *top and bottom* walls is *not* stored (those are
$y$-boundaries and $u$ lives at cell-centre height), so no-slip there **does**
need ghost cells: `u_ghost = 2*u_wall - u_adjacent`. Half your boundaries are
easy and half need ghosts, and which is which differs between `u` and `v`. This
is the bookkeeping; `code/ch06/mac_grid.m` builds it once so later tutorials just
call it.

### Interpolation: when you need $u$ where $v$ lives

The convective term needs products like $uv$ at cell **corners**. Interpolate:

$$u_{corner(i,j)} = \tfrac{1}{2}\left(u_{i,j-1} + u_{i,j}\right), \qquad
v_{corner(i,j)} = \tfrac{1}{2}\left(v_{i-1,j} + v_{i,j}\right)$$

Each is a simple average of the two nearest values — second-order accurate
because the corner is exactly midway. Tutorial 7 does this in two lines.

---

## 6.4 Poiseuille flow: your first validated solution

The flow between parallel plates separated by $h$, driven by $-dp/dx = G$:

$$\mu\frac{d^2u}{dy^2} = -G, \qquad u(\pm h/2) = 0$$

Integrating twice:

$$\boxed{u(y) = \frac{G}{2\mu}\left(\frac{h^2}{4} - y^2\right)}$$

with

$$u_{max} = \frac{Gh^2}{8\mu}, \qquad \bar{u} = \frac{Gh^2}{12\mu}, \qquad \frac{u_{max}}{\bar{u}} = \frac{3}{2}$$

That ratio $3/2$ is a fingerprint. **If your solver produces a parabola whose
peak is 1.5× its mean, the profile is right.** (For a circular pipe it's 2; for
a rectangular duct it's between the two and depends on aspect ratio — §6.5.)

The flow rate per unit depth, $Q = \bar{u}h = Gh^3/(12\mu)$, is the **cubic law**.
$Q \propto h^3$ at fixed pressure gradient means channel height dominates
everything: a 10% fabrication error in channel depth is a 30% flow-rate error.
This one scaling explains most of the frustration in experimental microfluidics.

`code/ch06/poiseuille.m` solves this on the MAC grid and reports the error, the
$u_{max}/\bar{u}$ ratio, flow rate vs the cubic law, and wall shear stress
$\tau_w = \mu\,du/dy|_{wall} = Gh/2$.

### A trap worth walking into: this is *not* exact

It is tempting to expect round-off agreement. The interior stencil has **zero**
truncation error on a quadratic, and the exact solution is a quadratic — so
surely the answer is exact?

It isn't. Run it and you get ~$2\times10^{-4}$ relative error at $n_y=64$. The
culprit is the cell-centred Dirichlet ghost. With the wall at $y=-h/2$, cell 1
centred at $y_1 = -h/2 + \Delta y/2$, and the ghost at $y_0 = -h/2-\Delta y/2$,
write $u = A(h^2/4 - y^2)$ and compare:

$$u(y_0) = A\left(-\tfrac{h\Delta y}{2} - \tfrac{\Delta y^2}{4}\right)
\qquad\text{but}\qquad
-u_1 = A\left(-\tfrac{h\Delta y}{2} + \tfrac{\Delta y^2}{4}\right)$$

The ghost relation $u_0 = -u_1$ is wrong by $A\Delta y^2/2$. That $O(\Delta y^2)$
boundary error is what you measure. The scheme is **second order, not exact**.

This generalizes, and it is the reason the exercise is here: *"my scheme is
second order" is a claim about the interior and the boundaries together, and the
boundary is almost always the weaker of the two.* So the right check is not "is
the error at round-off" but **"is the order 2"** — which is why `poiseuille.m`
now runs a grid-convergence study and reports the slope. Only claim exactness
when you have derived that the boundary treatment is exact too.

---

## 6.5 The 2D lie, and how big it is

Your simulations are 2D: infinitely deep channels. Real microchannels are
rectangular, typically wider than they are deep. How wrong is 2D?

For a rectangular duct of width $w$ and height $h$, the exact series solution
gives a mean velocity

$$\bar{u} = \frac{G h^2}{12\mu}\left[1 - \sum_{n,\text{odd}}\frac{192}{n^5\pi^5}\frac{h}{w}\tanh\frac{n\pi w}{2h}\right]$$

The bracket is the correction factor relative to the 2D result:

| $w/h$ | correction | 2D error |
|---|---|---|
| 1 (square) | 0.4217 | **137% too high** |
| 2 | 0.6860 | 46% too high |
| 5 | 0.8740 | 14% too high |
| 10 | 0.9370 | 7% too high |
| 20 | 0.9685 | 3% too high |
| ∞ | 1.0000 | exact |

**Read this table before you trust a 2D flow-rate prediction.** For a square
channel, 2D overpredicts flow rate by more than a factor of two — the two side
walls you deleted were doing half the work. 2D is a good model for a *shallow,
wide* channel ($w/h \gtrsim 10$) and a bad one for anything approaching square.

The saving grace: **2D is much better for flow *structure* than for flow
*magnitude*.** Recirculation patterns, streamline topology, where a droplet
pinches — these are mostly captured. Absolute flow rate and pressure drop are
not. So use 2D to understand behaviour, and correct the numbers with the table
above (or the empirical $f\!\cdot\!Re$ correlations) before quoting them.

`poiseuille.m` includes a `duct_correction(w/h)` helper that evaluates the series,
so you can apply this correction to any 2D result in the rest of the course.

---

## 6.6 Solving Stokes flow directly (and why we won't)

Because Stokes flow is linear, you *could* assemble one big matrix for
$(u, v, p)$ together and solve in a single shot:

$$\begin{pmatrix} \mu L & 0 & -G_x \\ 0 & \mu L & -G_y \\ D_x & D_y & 0\end{pmatrix}
\begin{pmatrix}u\\v\\p\end{pmatrix} = \begin{pmatrix}f_x\\f_y\\0\end{pmatrix}$$

This **monolithic** approach is exact, needs no timestepping, and is genuinely
attractive for pure Stokes problems. Two reasons this course goes a different
way:

1. The matrix is a **saddle-point system** — indefinite, with zeros on the
   diagonal. `\` handles it at these sizes, but it needs specialized
   preconditioners to scale, and it is not a gentle introduction.
2. It only works when the equations are linear. Add inertia (finite $Re$), a
   moving interface (Tutorial 11), or any nonlinearity, and it stops applying.

The **projection method** of Tutorial 7 handles the nonlinear case, reduces to a
sequence of ordinary Poisson solves you already know how to do, and is what
essentially every research incompressible-flow code uses. Learn that one.

---

## Exercises

**6.1** Run `checkerboard_demo`. Confirm the central-difference gradient of the
checkerboard is zero to machine precision, and that the staggered gradient is
not. Compute the ratio.

**6.2** Run `poiseuille`. Confirm the convergence study gives order 2 and that
$u_{max}/\bar{u} \to 1.5$ as the grid refines. Then derive a **second-order-exact**
boundary treatment: instead of the linear ghost $u_0 = -u_1$, fit a quadratic
through the wall value and the first two cell centres to get
$u_0 = -2u_1 + \tfrac{1}{3}u_2$ ... (work out the correct coefficients yourself).
Implement it and check whether the parabola now *is* reproduced to round-off.
*(It should be. This is the payoff for understanding where the error came from
rather than accepting "second order" as a brand.)*

**6.3** Verify the cubic law: compute $Q$ for $h$ = 25, 50, 100, 200 µm at fixed
$G$, and confirm the slope of $\log Q$ vs $\log h$ is 3.

**6.4** A design requires 10 µL/min through a 20 mm long, 200 µm wide channel.
Using the 2D result plus `duct_correction`, find the depth needed to keep the
pressure drop under 100 mbar. Then recompute with pure 2D and quantify how badly
you would have been misled.

**6.5** Confirm the $u_{max}/\bar{u}$ fingerprint distinguishes geometries:
compute it for your 2D solution (1.5), and evaluate the duct series for
$w/h = 1$. *(You should get ≈2.1 for the square duct centreline — the shape of
the profile changes, not just its magnitude.)*

**6.6 (grid-sense)** Print `size(u)`, `size(v)`, `size(p)` for `nx=8, ny=6`.
Write down, without running anything, which array index corresponds to the
physical wall at $x=0$, $x=L_x$, $y=0$, $y=L_y$ for each variable, and which of
those need ghost cells. Check against `mac_grid.m`. *(Do this. Every bug in
Tutorial 7 is a violation of this table.)*

---

**Next:** [Tutorial 7 — The projection method](07-projection-method.md).
Everything is now in place: staggered storage, a Poisson solver, advection
schemes, and stability limits. Time to write a real Navier–Stokes solver.
