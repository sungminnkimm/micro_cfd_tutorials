# Tutorial 4 — Unsteady Diffusion and Stability

**Goal:** Add time. Learn why your timestep is not a free parameter, derive the
stability limit rather than discovering it by crashing, and understand the
explicit/implicit trade-off that shapes every solver you write from here on.

**Code:** `code/ch04/diffusion1d.m`, `code/ch04/stability_demo.m`, `code/ch04/diffusion2d.m`

---

## 4.1 The equation

$$\frac{\partial c}{\partial t} = D\nabla^2 c$$

Physically: a dye spreading in still water, heat conducting, or — with $D\to\nu$
— momentum diffusing. This last reading matters: **the viscous term of
Navier–Stokes is exactly this equation**, so the stability limit you derive here
is the limit that will constrain your Tutorial 7 solver.

The exact solution for an initial point source of mass $M$ in 1D:

$$c(x,t) = \frac{M}{\sqrt{4\pi D t}}\exp\!\left(-\frac{x^2}{4Dt}\right)$$

a Gaussian whose width grows as $\sigma = \sqrt{2Dt}$. That $\sqrt{t}$ — not $t$
— is the mathematical statement of "diffusion is slow over long distances," and
it is why microfluidic mixing is hard. To go twice as far takes four times as
long.

---

## 4.2 Explicit (forward Euler): simple, and limited

Discretize time forward, space centrally:

$$\frac{c_i^{n+1} - c_i^n}{\Delta t} = D\frac{c_{i-1}^n - 2c_i^n + c_{i+1}^n}{\Delta x^2}$$

Solve for the new value — everything on the right is *known*:

$$c_i^{n+1} = c_i^n + r\left(c_{i-1}^n - 2c_i^n + c_{i+1}^n\right), \qquad \boxed{r = \frac{D\Delta t}{\Delta x^2}}$$

$r$ is the **diffusion number**. In MATLAB the whole timestep is one line:

```matlab
c(2:end-1) = c(2:end-1) + r*(c(1:end-2) - 2*c(2:end-1) + c(3:end));
```

or, with the matrix from Tutorial 2, `c = c + dt*D*(D2*c)`.

### Deriving the stability limit (von Neumann analysis)

This is worth doing once by hand, because the same argument gives you every
stability limit in the course.

Assume the error at step $n$ is a Fourier mode, $c_j^n = G^n e^{ikj\Delta x}$,
where $G$ is the per-step **amplification factor**. Stability requires
$|G| \le 1$ for every wavenumber $k$ — otherwise some mode grows without bound.

Substitute into the update:

$$G^{n+1}e^{ikj\Delta x} = G^n e^{ikj\Delta x}\left[1 + r\left(e^{-ik\Delta x} - 2 + e^{ik\Delta x}\right)\right]$$

Use $e^{i\theta} + e^{-i\theta} = 2\cos\theta$ and $1-\cos\theta = 2\sin^2(\theta/2)$:

$$G = 1 + r(2\cos k\Delta x - 2) = 1 - 4r\sin^2\!\left(\frac{k\Delta x}{2}\right)$$

$G$ is real and decreasing in $r$. It is worst for the highest wavenumber, where
$\sin^2 = 1$ — that is the **sawtooth mode** $(+1,-1,+1,-1,\dots)$, the one that
oscillates cell to cell. There $G = 1 - 4r$, and $|G|\le 1$ requires

$$\boxed{r = \frac{D\Delta t}{\Delta x^2} \le \frac{1}{2}}$$

In 2D the same argument with modes in both directions gives

$$D\,\Delta t\left(\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2}\right) \le \frac{1}{2}
\quad\Longrightarrow\quad
\Delta t \le \frac{1}{2D}\left(\frac{1}{\Delta x^2}+\frac{1}{\Delta y^2}\right)^{-1}$$

which on a uniform grid is $\Delta t \le \Delta x^2/(4D)$ — **half** the 1D limit.
Going to 3D costs you another third. Forgetting to tighten $\Delta t$ when moving
from 1D to 2D is a rite of passage; the symptom is instant, spectacular blowup.

### What instability looks like

At $r$ slightly over $1/2$, the sawtooth mode grows by $|1-4r|$ per step. At
$r = 0.6$ that's $1.4\times$ per step — after 100 steps, $10^{14}$. You will see:

1. A cell-to-cell zigzag appears, usually near a boundary or a sharp gradient.
2. It grows geometrically over tens of steps.
3. `NaN` or `Inf` everywhere.

**The zigzag is diagnostic.** A smooth-but-wrong answer is a *boundary condition*
bug. A cell-to-cell sawtooth that grows is a *stability* violation. Learning to
tell these apart on sight is worth more than any amount of extra reading.

### The cost, in microfluidic terms

Suppose you want to simulate mixing in a 100 µm channel with $D = 5\times10^{-10}$
m²/s, resolved with $\Delta x = 1$ µm:

$$\Delta t \le \frac{\Delta x^2}{4D} = \frac{10^{-12}}{2\times10^{-9}} = 5\times10^{-4}\ \text{s}$$

To reach the diffusive mixing time $L^2/D = 20$ s you need **40,000 steps**. Halve
$\Delta x$ and you need 160,000 (four times as many steps, each four times more
expensive in 2D — a 16× cost increase for a 2× resolution increase). This
$\Delta t \propto \Delta x^2$ scaling is what drives everyone to implicit methods.

---

## 4.3 Implicit (backward Euler): unconditionally stable

Evaluate the spatial derivative at the *new* time level:

$$\frac{c_i^{n+1} - c_i^n}{\Delta t} = D\frac{c_{i-1}^{n+1} - 2c_i^{n+1} + c_{i+1}^{n+1}}{\Delta x^2}$$

Now the unknowns are coupled — you must solve a linear system:

$$(I - D\Delta t\, L)\,c^{n+1} = c^n$$

where $L$ is the Laplacian matrix from Tutorial 2. In MATLAB:

```matlab
[~, L] = diff_matrix(N, dx, 'neumann0');
A = speye(N) - dt*D*L;
[Lf,Uf,Pf,Qf,Rf] = lu(A);              % once, outside the loop
...
c = Qf*(Uf\(Lf\(Pf*(Rf\c))));          % each step
```

Von Neumann analysis gives

$$G = \frac{1}{1 + 4r\sin^2(k\Delta x/2)}$$

which satisfies $0 < G \le 1$ for **every** $r > 0$. **Unconditionally stable.**
Take any timestep you like.

### But stable ≠ accurate

This is the point people miss. Backward Euler is first-order accurate in time,
$O(\Delta t)$. At $r = 100$ it will not blow up — it will just give you a heavily
damped, wrong answer, *smoothly and confidently*. Stability buys you the freedom
to choose $\Delta t$ on **accuracy** grounds instead of stability grounds. It does
not excuse you from choosing.

Practical rule: pick $\Delta t$ so the solution changes by a few percent per step,
then verify by halving it and checking the answer doesn't move.

---

## 4.4 Crank–Nicolson: second order in time

Average the two:

$$\left(I - \frac{D\Delta t}{2}L\right)c^{n+1} = \left(I + \frac{D\Delta t}{2}L\right)c^n$$

$$G = \frac{1 - 2r\sin^2(k\Delta x/2)}{1 + 2r\sin^2(k\Delta x/2)}$$

$|G| \le 1$ always — unconditionally stable — *and* second-order accurate in
time. Best of both, at the cost of one extra matrix-vector product.

**The caveat:** for large $r$, $G \to -1$ for the highest modes. Stable (magnitude
≤ 1) but *oscillatory* — sharp initial conditions produce ringing that decays
only slowly. If your initial condition has a discontinuity, either take a few
backward-Euler steps first to smooth it, or keep $r \lesssim 1$.

### Summary table

| Scheme | Stability | Time order | Cost/step | Use when |
|---|---|---|---|---|
| Forward Euler | $r \le 1/2$ | 1 | 1 matvec | $\Delta t$ already small for other reasons |
| Backward Euler | unconditional | 1 | 1 solve | steady-state hunting, stiff problems, rough ICs |
| Crank–Nicolson | unconditional | 2 | 1 solve + 1 matvec | accurate transients (**the default**) |
| RK3 / Adams–Bashforth | $r \lesssim 0.5$ | 3 / 2 | 3 / 1 matvec | explicit but need time accuracy (T7) |

---

## 4.5 The comparison script

`code/ch04/stability_demo.m` runs all three at $r = 0.4$ (stable), $r = 0.51$
(marginal), and $r = 0.6$ (unstable) so you can watch the sawtooth appear.

```matlab
%STABILITY_DEMO  Watch the explicit scheme cross its stability limit.
if ~exist('QUIET','var'), QUIET = false; end

L = 100e-6;  D = 5e-10;  N = 101;
dx = L/(N-1);  x = linspace(0,L,N)';
c0 = double(abs(x - L/2) < L/10);        % a top-hat of dye

rvals  = [0.40 0.51 0.60];
nsteps = 400;

for q = 1:numel(rvals)
    r  = rvals(q);
    dt = r*dx^2/D;
    c  = c0;
    for n = 1:nsteps
        c(2:end-1) = c(2:end-1) + r*(c(1:end-2) - 2*c(2:end-1) + c(3:end));
        c(1) = c(2);  c(end) = c(end-1);      % zero-flux walls
        if ~all(isfinite(c)), break; end
    end
    % A sawtooth mode has large cell-to-cell curvature: measure it directly.
    % diff(c,2) is the undivided second difference -- big for (+1,-1,+1,...).
    if all(isfinite(c)), saw = max(abs(diff(c,2))); else, saw = Inf; end
    if r <= 0.5, tag = 'STABLE'; else, tag = 'UNSTABLE'; end
    fprintf('  r = %.2f : max|c| = %10.3e   sawtooth = %10.3e   %s\n', ...
            r, max(abs(c)), saw, tag);
end
```

Run it. At $r = 0.40$ the top hat spreads into a smooth Gaussian. At $r = 0.51$
— *one percent* over the limit — it survives 400 steps but the sawtooth
indicator has grown measurably. At $r = 0.60$ everything is `Inf`.

That $r=0.51$ case is the important one. **Stability limits are sharp, not
soft.** There is no "slightly unstable"; there is only "growing slowly enough
that you haven't noticed yet."

---

## 4.6 Validating against the exact solution

`code/ch04/diffusion1d.m` solves with all three schemes and compares to the
analytic Gaussian, giving you a **time-order** convergence study to go with the
space-order study from Tutorial 2.

The procedure mirrors Tutorial 2's, with one wrinkle: to isolate the *time*
error you must hold the grid fixed and refine only $\Delta t$. If you refine both
at once you measure some mixture of the two and get a meaningless number.

```matlab
% Fixed fine grid so spatial error is negligible; refine dt only.
dts = dt0 ./ 2.^(0:4);
for k = 1:numel(dts)
    c = march(c0, dts(k), scheme);
    err(k) = norm(c - c_exact(tEnd), inf);
end
p = log2(err(1:end-1)./err(2:end));       % expect 1, 1, 2
```

Expected: forward Euler 1, backward Euler 1, Crank–Nicolson 2.

**If Crank–Nicolson measures 1 instead of 2**, the usual cause is boundary
conditions applied at the wrong time level — CN needs them averaged between $n$
and $n+1$ too. Applying BCs only at $n+1$ silently drops you to first order. This
is a real and very common bug, and the convergence study is the only thing that
catches it.

---

## 4.7 2D, and what it costs

`code/ch04/diffusion2d.m` extends to 2D with the Tutorial 3 machinery. The
implicit operator is $(I - D\Delta t\,L_{2D})$, and — the key structural point —
you factorize it **once**:

```matlab
Lap = kron(speye(ny),Dxx) + kron(Dyy,speye(nx));
A   = speye(nx*ny) - dt*D*Lap;
[Lf,Uf,Pf,Qf,Rf] = lu(A);          % ONCE
for n = 1:nsteps
    c = Qf*(Uf\(Lf\(Pf*(Rf\c(:)))));
    c = reshape(c,nx,ny);
end
```

This only works because $\Delta t$ and $D$ are **constant**. If you vary the
timestep adaptively, you must refactorize — which usually costs more than it
saves. Practical advice: pick a fixed $\Delta t$ and keep it.

### The operator-splitting alternative (ADI)

For large 2D/3D grids, factorizing the full 2D operator is expensive. **Alternating
Direction Implicit** splits the step into two half-steps, each implicit in only
one direction:

$$\left(I - \tfrac{D\Delta t}{2}L_x\right)c^* = \left(I + \tfrac{D\Delta t}{2}L_y\right)c^n$$
$$\left(I - \tfrac{D\Delta t}{2}L_y\right)c^{n+1} = \left(I + \tfrac{D\Delta t}{2}L_x\right)c^*$$

Each half-step is a set of independent **tridiagonal** solves, which cost $O(N)$
instead of $O(N^{1.5})$. Unconditionally stable, second-order. Worth knowing
about; not needed at the grid sizes in this course.

---

## 4.8 What carries into Navier–Stokes

Tutorial 7's solver will treat the viscous term explicitly, so this limit binds:

$$\Delta t \le \frac{1}{2\nu}\left(\frac{1}{\Delta x^2}+\frac{1}{\Delta y^2}\right)^{-1}$$

with $\nu = 10^{-6}$ m²/s for water — **2000× larger than the mass diffusivity
$D$**. That factor is the Schmidt number from Tutorial 1, and here is what it
means operationally:

> The viscous stability limit is 2000× more restrictive than the scalar
> diffusion limit. Momentum is what constrains your timestep, not the dye you
> actually care about watching.

Concretely, at $\Delta x = 1$ µm in water: $\Delta t \le 2.5\times10^{-7}$ s. To
simulate 1 second of mixing, that's 4 million steps. This is *the* reason
Tutorial 9 solves the steady flow field once, freezes it, and then transports the
scalar with its own much larger timestep. Understanding why that shortcut is
legitimate — and it is, because $Re \ll 1$ makes the flow steady — is one of the
most valuable things in this course.

---

## Exercises

**4.1** Run `stability_demo`. Confirm the $r=0.51$ case is measurably worse than
$r=0.40$ despite both "running". Now run $r=0.51$ for 5000 steps. Does it survive?

**4.2** Derive the amplification factor for the **2D** explicit scheme by
assuming $c^n_{jk} = G^n e^{i(k_x j\Delta x + k_y k\Delta y)}$. Confirm the factor
of 1/4 rather than 1/2 on a uniform grid.

**4.3** Run the time-convergence study in `diffusion1d`. Confirm orders 1, 1, 2.
Then deliberately apply the Crank–Nicolson boundary condition only at level
$n+1$ and re-measure. Watch it drop to 1.

**4.4** Set up a 100 µm channel cross-section with dye at concentration 1 in the
left half, 0 in the right, zero-flux walls. How long until the concentration is
uniform to within 1%? Compare to the $L^2/D$ estimate from Tutorial 1.

**4.5 (the microfluidics one)** Two streams meet in a Y-junction and flow side by
side down a 50 µm channel at 2 mm/s. **Model this as a 1D diffusion problem in
the cross-stream direction, where time maps to downstream distance via
$x = Ut$.** How far down the channel before they're mixed? Compare to
$L_{mix} = Pe\cdot w$ from Tutorial 1. *(This trick — trading a steady 2D problem
for an unsteady 1D one — is exact when diffusion along the flow is negligible,
i.e. at high Pe. You'll check it against a real 2D solve in Tutorial 9.)*

**4.6** Time the explicit and implicit 2D solvers over the same physical
duration, on grids of 32², 64², 128², 256². Explicit needs $\Delta t \propto \Delta x^2$;
implicit can hold $\Delta t$ fixed. Plot total runtime vs $N$ for both and find
the crossover grid size. *(This is the measurement that should decide the
question for your own problem — not folklore.)*

---

**Next:** [Tutorial 5 — Advection and the Péclet number](05-advection-and-peclet.md).
Diffusion is well-behaved and forgiving. Advection is neither, and at $Pe = 200$
advection is what you have.
