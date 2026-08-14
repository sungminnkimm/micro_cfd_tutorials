# Tutorial 2 — Finite Differences and Truncation Error

**Goal:** Turn derivatives into matrix operations on a grid, and learn the single
most important debugging technique in CFD: the **grid convergence study**.

**Code:** `code/ch02/fd_stencils.m`, `code/ch02/convergence_study.m`, `code/ch02/diff_matrix.m`

---

## 2.1 The whole idea in one line

A CFD code is a machine that replaces $\frac{\partial}{\partial x}$ with a
weighted sum of neighbouring grid values. That's it. Everything else — grids,
solvers, boundary conditions — is plumbing around that substitution.

The substitution is *approximate*. The error it makes is called **truncation
error**, it depends on the grid spacing $\Delta x$, and if you cannot predict
that dependence you cannot tell a converged solution from a wrong one.

---

## 2.2 Deriving a stencil from Taylor series

Take a smooth $f(x)$ and expand about $x_i$ with spacing $h = \Delta x$:

$$f_{i+1} = f_i + h f'_i + \frac{h^2}{2}f''_i + \frac{h^3}{6}f'''_i + O(h^4)$$
$$f_{i-1} = f_i - h f'_i + \frac{h^2}{2}f''_i - \frac{h^3}{6}f'''_i + O(h^4)$$

**Forward difference.** Rearrange the first:

$$f'_i = \frac{f_{i+1} - f_i}{h} - \underbrace{\frac{h}{2}f''_i + O(h^2)}_{\text{truncation error}}$$

Leading error term $\propto h^1$. **First-order accurate.**

**Central difference.** Subtract the second from the first — the $f''$ terms
cancel:

$$f'_i = \frac{f_{i+1} - f_{i-1}}{2h} - \frac{h^2}{6}f'''_i + O(h^4)$$

Leading error $\propto h^2$. **Second-order accurate.** Same cost as forward
differencing, one order better. This is why central differencing is the default.

**Second derivative.** Add them instead — the $f'$ and $f'''$ terms cancel:

$$f''_i = \frac{f_{i+1} - 2f_i + f_{i-1}}{h^2} - \frac{h^2}{12}f''''_i + O(h^4)$$

Second-order accurate. This three-point stencil `[1 -2 1]/h²` is the workhorse of
this entire course — it is the discrete Laplacian, and it appears in the
diffusion term, the pressure Poisson equation, and the viscous term of
Navier–Stokes.

### What "second-order" actually promises

If a scheme is $p$-th order, halving $h$ divides the error by $2^p$:

| Order | Halve the grid → error drops by | Need 10× less error → refine by |
|-------|--------------------------------|--------------------------------|
| 1st | 2× | 10× (100× the cells in 2D) |
| 2nd | 4× | 3.2× (10× the cells in 2D) |
| 4th | 16× | 1.8× (3.2× the cells in 2D) |

Second order is the sweet spot for CFD and it is what everything in this course
delivers. **If you measure something other than 2, you have a bug.** That claim
is the reason this tutorial exists.

---

## 2.3 Measuring the order — the convergence study

The procedure, which you will repeat in every remaining tutorial:

1. Pick a problem with a known exact answer.
2. Solve it on grids of spacing $h, h/2, h/4, h/8, \dots$
3. Compute the error norm $E(h) = \|f_{\text{numerical}} - f_{\text{exact}}\|$.
4. Plot $\log E$ vs $\log h$. **The slope is the observed order of accuracy.**
5. Compare to the theoretical order. They must match.

Step 5 is the assertion. If theory says 2 and you measure 1, something is
first-order-wrong — usually a boundary condition, occasionally an index slip that
shifts the stencil by half a cell.

Here is the study, `code/ch02/convergence_study.m`:

```matlab
%CONVERGENCE_STUDY  Measure the observed order of accuracy of FD stencils.
if ~exist('QUIET','var'), QUIET = false; end

% Test function on [0,1] with known derivatives.  Choose something with all
% derivatives nonzero so no error term accidentally vanishes.
f   = @(x) sin(4*x) .* exp(-x);
df  = @(x) 4*cos(4*x).*exp(-x) - sin(4*x).*exp(-x);
d2f = @(x) -16*sin(4*x).*exp(-x) - 8*cos(4*x).*exp(-x) + sin(4*x).*exp(-x);

N = 2.^(4:11);                 % 16 .. 2048 intervals
h = 1./N;
[e_fwd, e_cen, e_2nd] = deal(zeros(size(N)));

for k = 1:numel(N)
    x = linspace(0,1,N(k)+1)';
    u = f(x);

    % --- interior nodes only: boundaries need one-sided stencils, and mixing
    % --- them in would pollute the order we are trying to measure.
    i = 2:N(k);

    d_fwd = (u(i+1) - u(i)) / h(k);                    % O(h)
    d_cen = (u(i+1) - u(i-1)) / (2*h(k));              % O(h^2)
    d_2nd = (u(i+1) - 2*u(i) + u(i-1)) / h(k)^2;       % O(h^2)

    % Max-norm error.  (Use norm(...,inf); the 2-norm needs a sqrt(h) weight
    % to be a true function norm -- a classic way to fake a good-looking rate.)
    e_fwd(k) = norm(d_fwd - df(x(i)),  inf);
    e_cen(k) = norm(d_cen - df(x(i)),  inf);
    e_2nd(k) = norm(d_2nd - d2f(x(i)), inf);
end

% Observed order between successive grids: p = log2(E(h)/E(h/2))
p_fwd = log2(e_fwd(1:end-1)./e_fwd(2:end));
p_cen = log2(e_cen(1:end-1)./e_cen(2:end));
p_2nd = log2(e_2nd(1:end-1)./e_2nd(2:end));

fprintf('\n     N   fwd err   order |  cen err   order |  2nd err   order\n');
fprintf('  ------------------------------------------------------------\n');
for k = 1:numel(N)
    fprintf('  %4d  %9.2e', N(k), e_fwd(k));
    if k>1, fprintf(' %6.2f |', p_fwd(k-1)); else, fprintf('    --- |'); end
    fprintf(' %9.2e', e_cen(k));
    if k>1, fprintf(' %6.2f |', p_cen(k-1)); else, fprintf('    --- |'); end
    fprintf(' %9.2e', e_2nd(k));
    if k>1, fprintf(' %6.2f\n', p_2nd(k-1)); else, fprintf('    ---\n'); end
end
fprintf('\n  Asymptotic orders: fwd %.2f | cen %.2f | 2nd %.2f\n', ...
        p_fwd(end), p_cen(end), p_2nd(end));

if ~QUIET
    figure('Color','w','Position',[100 100 620 480]);
    loglog(h, e_fwd,'o-','LineWidth',1.6); hold on
    loglog(h, e_cen,'s-','LineWidth',1.6);
    loglog(h, e_2nd,'^-','LineWidth',1.6);
    loglog(h, h,      'k--');
    loglog(h, h.^2, 'k:');
    grid on; xlabel('h = \Deltax'); ylabel('max-norm error');
    legend({'forward d/dx','central d/dx','central d^2/dx^2','O(h)','O(h^2)'}, ...
           'Location','southeast');
    title('Grid convergence: slope = observed order');
end
```

Run it. You should see the forward-difference column converge to 1.00 and both
central columns to 2.00.

### Two traps this script deliberately avoids

**Trap 1 — including the boundary.** If you naively apply the central stencil at
$i=1$ you index outside the array. If you patch it with a first-order one-sided
formula, your *global* error becomes first order and the whole study reports 1
instead of 2. That is an extremely common and very confusing bug. Here we measure
interior nodes only, and revisit boundaries properly in Tutorial 3.

**Trap 2 — the wrong norm.** `norm(e)` in MATLAB is the vector 2-norm,
$\sqrt{\sum e_i^2}$. As you refine, the number of points grows, so this
*mechanically* grows even at constant pointwise error. The function-space
$L_2$ norm is $\sqrt{h\sum e_i^2}$ — note the $\sqrt{h}$. Either use `inf`
(max) norm as here, or remember the weight. Reporting an unweighted 2-norm is a
reliable way to convince yourself a broken code is second order.

### When the rate is right but not 2

- **Rate = 1 everywhere:** boundary treatment, or a one-sided stencil leaking in.
- **Rate = 2 then drops toward 0 on the finest grids:** you have hit round-off.
  At double precision the error floor is $\sim \varepsilon/h^2 \approx 10^{-16}/h^2$
  for a second derivative. At $h=10^{-5}$ that floor is $10^{-6}$. This is real
  and expected — stop refining.
- **Rate ≈ 1.5, non-integer:** your exact solution isn't smooth enough, or you
  are measuring in a region containing a kink/discontinuity.
- **Rate = 4 when you expected 2:** you got lucky with a symmetric test function
  whose $f'''$ vanishes. Change the test function.

---

## 2.4 Derivatives as matrices

Writing `(u(i+1)-u(i-1))/(2*h)` inside a loop is fine for one derivative. But
implicit methods (Tutorial 4) and Poisson solves (Tutorial 3) need the derivative
as a **matrix** so you can invert it. Build it once with `spdiags`:

`code/ch02/diff_matrix.m`:

```matlab
function [D1, D2] = diff_matrix(N, h, bc)
%DIFF_MATRIX  Sparse 1-D first and second derivative operators.
%
%   [D1,D2] = diff_matrix(N,h,bc)  returns N-by-N sparse matrices acting on a
%   column vector of N nodal values with spacing h.
%
%   bc = 'dirichlet0'  u = 0 at both walls      -> ghost u_0 = -u_1
%      = 'neumann0'    du/dn = 0 at both walls  -> ghost u_0 = +u_1
%      = 'periodic'    wraps around

e  = ones(N,1);
D1 = spdiags([-e  0*e  e], [-1 0 1], N, N) / (2*h);
D2 = spdiags([ e -2*e  e], [-1 0 1], N, N) / h^2;

switch lower(bc)
    case 'periodic'
        D1(1,N) = -1/(2*h);   D1(N,1) =  1/(2*h);
        D2(1,N) =  1/h^2;     D2(N,1) =  1/h^2;

    case 'dirichlet0'
        D2(1,1) = -3/h^2;     D2(N,N) = -3/h^2;
        D1(1,1) =  1/(2*h);   D1(N,N) = -1/(2*h);

    case 'neumann0'
        D2(1,1) = -1/h^2;     D2(N,N) = -1/h^2;
        D1(1,1) = -1/(2*h);   D1(N,N) =  1/(2*h);

    otherwise
        error('diff_matrix:bc','Unknown bc ''%s''.', bc);
end
end
```

**Derive those corner entries yourself — do not trust mine.** Cell 1 is centred
at $x = h/2$, so the wall at $x=0$ sits half a cell to its left, where a ghost
value $u_0$ lives. For a Dirichlet wall at $u=0$, linear interpolation between
$u_0$ and $u_1$ must vanish at the wall, so $u_0 = -u_1$. Substitute into the
standard stencils:

$$\left.\frac{\partial^2 u}{\partial x^2}\right|_1 = \frac{u_0 - 2u_1 + u_2}{h^2} = \frac{-3u_1 + u_2}{h^2}
\qquad
\left.\frac{\partial u}{\partial x}\right|_1 = \frac{u_2 - u_0}{2h} = \frac{+u_1 + u_2}{2h}$$

which is where `-3/h^2` and `+1/(2h)` come from. For Neumann, $u_0 = +u_1$ and
the signs flip. **The sign on the `D1` corner entry is the single easiest thing
in this file to get backwards** — it flips between the two boundary types and
between the left and right walls, and the resulting bug is subtle (the solution
looks plausible, the convergence rate quietly drops to 1). Exercise 2.6 catches it.

One accuracy caveat: with these ghost relations `D2` is globally **second order**
— that is the standard finite-volume result and the payoff for cell-centring.
The boundary *rows* of `D1` are only first order. This course never
differentiates a scalar right at a wall, so it doesn't bite, but don't reuse
those corner rows in a high-accuracy setting.

Now `D2*u` *is* $\partial^2 u/\partial x^2$, and — critically — `(I - dt*nu*D2)\u`
is an implicit diffusion step. Everything downstream depends on being able to
write derivatives this way.

**Always look at your matrix before trusting it:**

```matlab
[D1,D2] = diff_matrix(8, 0.1, 'dirichlet0');
full(D2)          % eyeball the numbers
spy(D2)           % eyeball the sparsity pattern
```

The three-band structure should be obvious, and the corners should show your
boundary modification. Ninety percent of "my solver diverges" turns out to be a
matrix that doesn't look like this.

---

## 2.5 Node-centred vs cell-centred: pick one and commit

Two conventions for laying $N$ unknowns on a domain $[0,L]$:

**Node-centred** — unknowns *on* the boundaries:
```
 x:  0     h    2h   ...        L
     |-----|-----|-----|-----|--|
     u1    u2    u3          uN
```
`x = linspace(0,L,N)`, `h = L/(N-1)`. A Dirichlet condition is just `u(1)=value`.
Natural for finite differences.

**Cell-centred** — unknowns at cell midpoints, boundaries fall *between* nodes:
```
 x:  0                          L
     |--x--|--x--|--x--|--x--|--x--|
       u1    u2    u3    u4    u5
```
`x = (h/2 : h : L-h/2)`, `h = L/N`. A Dirichlet condition needs a ghost value
`u_0 = 2*u_wall - u_1`. Natural for finite volumes and for **conservation**,
because each unknown owns a control volume whose fluxes you can balance exactly.

**This course uses cell-centred** from Tutorial 3 onward, because the staggered
(MAC) grid in Tutorial 6 is built on it and because mass conservation is
something you want to be exact, not approximate. Mixing the two conventions in
one code is a rich source of half-cell offset bugs — the symptom is a solution
that looks right but converges at order 1.

---

## 2.6 A note on stencils you'll meet later

| Stencil | Formula | Order | Where it appears |
|---|---|---|---|
| Central 1st | $(f_{i+1}-f_{i-1})/2h$ | 2 | pressure gradient, low-Re convection |
| Central 2nd | $(f_{i+1}-2f_i+f_{i-1})/h^2$ | 2 | all diffusion, all Poisson |
| Upwind 1st | $(f_i-f_{i-1})/h$ for $u>0$ | 1 | high-Pe scalar transport (T5) |
| Face interpolation | $f_{i+1/2}=(f_i+f_{i+1})/2$ | 2 | staggered grid (T6) |
| One-sided 2nd | $(-3f_i+4f_{i+1}-f_{i+2})/2h$ | 2 | outflow BCs (T8) |

You do not need to memorize these. You need to know that each has an order, and
that mixing a first-order one into a second-order code drags the whole thing to
first order.

---

## Exercises

**2.1** Run `convergence_study`. Confirm the orders are 1, 2, 2. Now change the
error norm from `inf` to plain `norm(...)` (unweighted 2-norm) and rerun. What
orders do you measure now? Explain the shift. *(This is Trap 2, and you should
see it happen with your own eyes once.)*

**2.2** Push `N` to `2.^(4:16)`. Plot the second-derivative error. Find the grid
size where round-off takes over and the error starts *rising*. Verify it's near
$h \approx \varepsilon^{1/4}$ ≈ $1.2\times10^{-4}$.

**2.3** Derive the fourth-order central first derivative,
$f'_i = \frac{f_{i-2} - 8f_{i-1} + 8f_{i+1} - f_{i+2}}{12h}$, by combining Taylor
expansions to kill the $h^2$ term. Add it to the convergence study and confirm
slope 4.

**2.4** Build `diff_matrix(64, h, 'periodic')` and compute the eigenvalues of
`D2`. Confirm they are all real and non-positive, and that the largest-magnitude
one is close to $-4/h^2$. *(You will need that number in Tutorial 4 — it sets the
explicit diffusion timestep.)*

**2.5** Write a `diff_matrix2d(nx, ny, hx, hy, bc)` that builds the 2D Laplacian
using `kron`:
```matlab
L2D = kron(speye(ny), D2x) + kron(D2y, speye(nx));
```
Verify it on $f = \sin(\pi x)\sin(\pi y)$, whose Laplacian is $-2\pi^2 f$.
Confirm second order. *(You'll use this directly in Tutorial 3 — do it now.)*

**2.6 (the sign check)** Verify the `D1` corner rows numerically instead of
trusting the derivation. Take $u(x) = x$ on $[0,1]$ sampled at cell centres.
With `'neumann0'`, `D1*u` should be $\approx 1$ in the interior — but check what
it gives in row 1, and confirm it matches the ghost relation $u_0 = u_1$ (i.e.
$(u_2-u_1)/2h = 1/2$, not 1). Then deliberately flip the sign of `D1(1,1)` and
rerun `convergence_study` with a scheme that uses it. Watch the measured order
fall from 2 to 1. **This is what a boundary bug looks like**, and recognizing the
signature — right answer in the middle, order-1 convergence overall — will save
you hours in Tutorial 7.

---

**Next:** [Tutorial 3 — Boundary conditions and sparse linear systems](03-poisson-and-sparse-systems.md).
You can now differentiate. Next you learn to *invert*, which is what solving a
PDE actually means.
