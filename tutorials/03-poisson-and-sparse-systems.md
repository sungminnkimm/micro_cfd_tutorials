# Tutorial 3 — Boundary Conditions and Sparse Linear Systems

**Goal:** Solve $\nabla^2 \phi = f$ on a 2D grid. This is not a detour — the
pressure Poisson equation *is* this equation, and it will consume ~80% of the
runtime of your Navier–Stokes solver. Get it right here, once.

**Code:** `code/ch03/poisson2d.m`, `code/ch03/poisson_demo.m`, `code/ch03/poisson_mms.m`

---

## 3.1 Why Poisson matters so much

Three of the four things you will build in this course are Poisson solves in
disguise:

| Problem | Equation |
|---|---|
| Steady heat conduction | $\nabla^2 T = -q/k$ |
| Pressure in the projection method (T7) | $\nabla^2 \phi = \frac{1}{\Delta t}\nabla\cdot\mathbf{u}^*$ |
| Stream function from vorticity | $\nabla^2\psi = -\omega$ |
| Implicit diffusion step (T4) | $(I - \Delta t\, \nu\nabla^2)u^{n+1} = u^n$ |

Learn one solver, reuse it four times.

---

## 3.2 The 2D discrete Laplacian

On a cell-centred grid with spacing $\Delta x, \Delta y$, at cell $(i,j)$:

$$
\nabla^2\phi \approx
\frac{\phi_{i-1,j} - 2\phi_{i,j} + \phi_{i+1,j}}{\Delta x^2}
+ \frac{\phi_{i,j-1} - 2\phi_{i,j} + \phi_{i,j+1}}{\Delta y^2}
$$

the **five-point stencil**:

```
                 1/dy²
                   |
      1/dx² ---  -2/dx² - 2/dy²  --- 1/dx²
                   |
                 1/dy²
```

To turn this into a matrix we must flatten the 2D array of unknowns into a
vector. **Choose a convention and never deviate.** This course uses MATLAB's
native column-major order:

$$k = i + (j-1)\,n_x, \qquad i = 1\dots n_x,\ j = 1\dots n_y$$

which is exactly what `reshape(PHI, [], 1)` and `PHI(:)` do. So:

- `k-1`, `k+1` are the $x$-neighbours (`i∓1`)
- `k-nx`, `k+nx` are the $y$-neighbours (`j∓1`)

and to get back, `PHI = reshape(phi, nx, ny)`. Because MATLAB is column-major
and we put $i$ first, **`PHI(i,j)` means `PHI(x-index, y-index)`** — so when you
plot, you will need `PHI.'` (transpose) to get $x$ horizontal. Every plotting
call in this course does that; if your results look transposed, this is why.

### Building it with `kron`

Assembling the matrix with explicit loops over `k` is error-prone. The
Kronecker-product construction is three lines and hard to get wrong:

$$L_{2D} = I_{n_y} \otimes D_{2x} \; + \; D_{2y} \otimes I_{n_x}$$

```matlab
Lap = kron(speye(ny), D2x) + kron(D2y, speye(nx));
```

Read it as: "the $x$-derivative acts within each row of cells and is repeated
$n_y$ times; the $y$-derivative connects cells $n_x$ apart." The ordering of the
two `kron` arguments encodes the column-major flattening. If you swap them you
get the transposed grid — a bug that produces a *plausible-looking but wrong*
answer on a square grid and an immediate size error on a rectangular one.
**Always test new grid code on a non-square grid** ($n_x \ne n_y$); it catches
this class of bug instantly.

---

## 3.3 Boundary conditions, properly

Cell-centred means no unknown sits on the boundary, so all BCs act through ghost
cells. Let the wall be between cell 1 and ghost cell 0.

### Dirichlet: $\phi = \phi_w$ at the wall

Linear interpolation to the wall must give $\phi_w$:

$$\frac{\phi_0 + \phi_1}{2} = \phi_w \quad\Rightarrow\quad \phi_0 = 2\phi_w - \phi_1$$

Substituting into the stencil at cell 1:

$$\frac{\phi_0 - 2\phi_1 + \phi_2}{\Delta x^2} = \frac{-3\phi_1 + \phi_2}{\Delta x^2} + \underbrace{\frac{2\phi_w}{\Delta x^2}}_{\text{goes to RHS}}$$

So Dirichlet does two things: it changes the diagonal from $-2$ to $-3$, **and it
adds a known term to the right-hand side.** Forgetting the RHS contribution is
the classic error — you get the homogeneous solution, which for $\phi_w=0$ looks
completely fine, so the bug hides until the day you use a nonzero wall value.

### Neumann: $\partial\phi/\partial n = g$ at the wall

$$\frac{\phi_1 - \phi_0}{\Delta x} = g \quad\Rightarrow\quad \phi_0 = \phi_1 - g\Delta x$$

$$\frac{\phi_0 - 2\phi_1 + \phi_2}{\Delta x^2} = \frac{-\phi_1 + \phi_2}{\Delta x^2} - \underbrace{\frac{g}{\Delta x}}_{\text{to RHS}}$$

Diagonal goes from $-2$ to $-1$.

### The all-Neumann trap — read this twice

If **every** boundary is Neumann, the matrix is **singular**. Physically: if you
only specify gradients, the solution is determined only up to an additive
constant — $\phi$ and $\phi + 42$ both satisfy every equation. Numerically:
`A` has a null space spanned by the constant vector `ones(N,1)`, so `A\b` will
either warn about singularity, return garbage, or return something huge.

**This is exactly the situation in the pressure Poisson equation for a closed
domain**, so you will meet it in Tutorial 7. Two fixes:

1. **Pin one value.** Replace one row with $\phi_1 = 0$. Simple, robust, and what
   this course does. It makes the matrix nonsingular at the cost of one
   equation's accuracy — irrelevant, because only $\nabla\phi$ is ever used.
2. **Project out the null space.** Subtract the mean from `b` and from the
   solution. Mathematically cleaner, required for iterative solvers.

There is also a **compatibility condition**: for a solution to exist at all,

$$\int_\Omega f\,dV = \oint_{\partial\Omega} g\,dS$$

("total source = total flux out"). If your right-hand side violates it, no
solution exists and the solver returns nonsense no matter how you regularize. In
the projection method this condition is equivalent to *net mass flux in = net
mass flux out*, which is why an inconsistent outflow BC makes the pressure solve
blow up. Remember that sentence; it will save you in Tutorial 8.

---

## 3.4 The solver

`code/ch03/poisson2d.m` builds and factorizes the operator. Note the split
between **assembly** (expensive, do once) and **solve** (cheap, do every
timestep) — this structure is what makes the Tutorial 7 solver fast.

```matlab
function S = poisson2d(nx, ny, dx, dy, bcTypes)
%POISSON2D  Assemble and factorize the 2-D cell-centred Laplacian.
%
%   S   = poisson2d(nx,ny,dx,dy,{'D','D','N','N'})   % {xlo,xhi,ylo,yhi}
%   phi = poisson_solve(S, f, bcVals)
%
%   'D' = Dirichlet, 'N' = Neumann.  bcVals is a 1x4 cell of wall values
%   (scalars or vectors along the wall).

ex = ones(nx,1);  ey = ones(ny,1);
Dxx = spdiags([ex -2*ex ex],[-1 0 1],nx,nx)/dx^2;
Dyy = spdiags([ey -2*ey ey],[-1 0 1],ny,ny)/dy^2;

% Boundary modifications to the diagonal.
%   Dirichlet: -2 -> -3     Neumann: -2 -> -1
adj = @(t) (t=='D')*(-1) + (t=='N')*(+1);
Dxx(1,1)   = Dxx(1,1)   + adj(bcTypes{1})/dx^2;
Dxx(nx,nx) = Dxx(nx,nx) + adj(bcTypes{2})/dx^2;
Dyy(1,1)   = Dyy(1,1)   + adj(bcTypes{3})/dy^2;
Dyy(ny,ny) = Dyy(ny,ny) + adj(bcTypes{4})/dy^2;

A = kron(speye(ny), Dxx) + kron(Dyy, speye(nx));

% All-Neumann -> singular.  Pin cell 1 to remove the constant null mode.
S.pinned = all(strcmp(bcTypes,'N'));
if S.pinned
    A(1,:) = 0;  A(1,1) = 1;
end

S.nx=nx; S.ny=ny; S.dx=dx; S.dy=dy; S.bcTypes=bcTypes;
S.A = A;
[S.L, S.U, S.P, S.Q, S.R] = lu(A);    % factorize ONCE
end
```

The solve applies the boundary right-hand-side contributions and back-substitutes
(`code/ch03/poisson_solve.m`):

```matlab
function phi = poisson_solve(S, f, bcVals)
b = f;                                        % nx-by-ny source array

% For a 'D' wall, bcVals is the wall VALUE phi_w.
% For an 'N' wall, bcVals is the OUTWARD normal derivative g = dphi/dn.
if S.bcTypes{1}=='D', b(1,:)   = b(1,:)   - 2*bcVals{1}(:).'/S.dx^2;
else,                 b(1,:)   = b(1,:)   -   bcVals{1}(:).'/S.dx;   end
if S.bcTypes{2}=='D', b(end,:) = b(end,:) - 2*bcVals{2}(:).'/S.dx^2;
else,                 b(end,:) = b(end,:) -   bcVals{2}(:).'/S.dx;   end
if S.bcTypes{3}=='D', b(:,1)   = b(:,1)   - 2*bcVals{3}(:)/S.dy^2;
else,                 b(:,1)   = b(:,1)   -   bcVals{3}(:)/S.dy;     end
if S.bcTypes{4}=='D', b(:,end) = b(:,end) - 2*bcVals{4}(:)/S.dy^2;
else,                 b(:,end) = b(:,end) -   bcVals{4}(:)/S.dy;     end

b = b(:);
if S.pinned, b(1) = 0; end                    % match the pinned row

phi = S.Q * (S.U \ (S.L \ (S.P * (S.R \ b))));
phi = reshape(phi, S.nx, S.ny);
end
```

### Every sign is a minus — and that is a *check*, not a coincidence

Derive the Neumann case for the left wall. The outward normal there is $-\hat{x}$,
so $g = \partial\phi/\partial n = -\partial\phi/\partial x$, giving
$\phi_0 = \phi_1 + g\Delta x$. Substitute:

$$\frac{\phi_0 - 2\phi_1 + \phi_2}{\Delta x^2} = \frac{-\phi_1+\phi_2}{\Delta x^2} + \frac{g}{\Delta x}$$

Move the known part right: $b_1 = f_1 - g/\Delta x$. Now do the **right** wall,
where the outward normal is $+\hat{x}$ and $\phi_{n_x+1} = \phi_{n_x} + g\Delta x$;
you get $b_{n_x} = f_{n_x} - g/\Delta x$ — *also* minus.

That symmetry is the point. Because $g$ is referred to the **outward** normal at
every wall, both walls contribute with the same sign. Sum every equation over the
whole grid: all interior flux terms cancel in pairs, and what survives is

$$\int_\Omega f\,dV - \oint_{\partial\Omega} g\,dS = 0$$

which is precisely the compatibility condition from §3.3. **So if you get one of
these four signs wrong, an all-Neumann problem silently violates its own
solvability condition** and the solver returns a plausible-looking field that
satisfies nothing. Deriving the signs from the outward normal — rather than
pattern-matching them — is what keeps this consistent. Exercise 3.6 makes the
failure visible.

### Why factorize separately

```matlab
A\b                        % ~O(n^1.5) each call, refactorizes EVERY time
[L,U,P,Q,R] = lu(A);       % once
x = Q*(U\(L\(P*(R\b))));   % ~O(n log n) each call
```

In Tutorial 7 you call the pressure solve 10,000+ times with the *same matrix*
and a different RHS each step. Factorizing once is a 10–50× speedup on the whole
simulation. This is the single highest-leverage optimization in the course.

For grids beyond ~500×500, direct factorization runs out of memory and you switch
to an iterative method (`pcg` with an incomplete-Cholesky preconditioner) or a
DCT-based fast Poisson solver. Appendix A covers both. Up to 256×256 — which
covers every simulation in these tutorials — direct LU is faster and simpler.

---

## 3.5 Verifying it: the Method of Manufactured Solutions

You cannot check a Poisson solver against intuition. You check it against a
*manufactured* exact solution. The trick:

1. **Pick the answer first.** Say $\phi_{ex}(x,y) = \sin(\pi x)\cos(2\pi y)$.
2. **Differentiate it analytically** to find the source that produces it:
   $f = \nabla^2\phi_{ex} = -(\pi^2 + 4\pi^2)\sin(\pi x)\cos(2\pi y)$.
3. **Feed that $f$ to your solver** with BCs read off $\phi_{ex}$.
4. **Compare.** The error should converge at second order.

The beauty is that step 2 works for *any* function you can differentiate, so you
can manufacture a test for any operator, any BC combination, any geometry. This
is the professional standard for code verification and it is the subject of
Tutorial 13. `code/ch03/poisson_mms.m`:

```matlab
%POISSON_MMS  Verify poisson2d by the method of manufactured solutions.
if ~exist('QUIET','var'), QUIET = false; end

phi_ex = @(x,y) sin(pi*x).*cos(2*pi*y);
src    = @(x,y) -(pi^2 + 4*pi^2)*sin(pi*x).*cos(2*pi*y);

N = [16 32 64 128 256];
err = zeros(size(N));

for k = 1:numel(N)
    % Deliberately non-square: catches kron-ordering bugs.
    nx = N(k);  ny = round(N(k)*3/4);
    dx = 1/nx;  dy = 1/ny;
    xc = (dx/2 : dx : 1-dx/2)';        % cell centres
    yc = (dy/2 : dy : 1-dy/2)';
    [X,Y] = ndgrid(xc,yc);             % ndgrid, NOT meshgrid -- see below

    S = poisson2d(nx,ny,dx,dy,{'D','D','D','D'});

    % Dirichlet values evaluated ON the walls (x=0,1 and y=0,1)
    bcVals = { phi_ex(0,yc),  phi_ex(1,yc), ...
               phi_ex(xc,0),  phi_ex(xc,1) };

    phi = S.solve(src(X,Y), bcVals);
    err(k) = norm(phi(:) - reshape(phi_ex(X,Y),[],1), inf);
end

p = log2(err(1:end-1)./err(2:end));
fprintf('\n    nx    max error   order\n');
fprintf('  ------------------------------\n');
for k=1:numel(N)
    if k==1, fprintf('  %4d   %9.3e     ---\n', N(k), err(k));
    else,    fprintf('  %4d   %9.3e    %5.2f\n', N(k), err(k), p(k-1)); end
end
fprintf('\n  Observed order: %.2f  (expected 2.00)\n\n', p(end));
```

### `ndgrid` vs `meshgrid` — a real bug source

- `meshgrid(x,y)` returns arrays indexed `(row, col) = (y, x)`.
- `ndgrid(x,y)` returns arrays indexed `(x, y)`.

Our flattening convention is `(i,j) = (x,y)`, which matches **`ndgrid`**. Using
`meshgrid` silently transposes your source term. On a square grid with a
symmetric test function it still "works", which is the worst possible failure
mode. **This course uses `ndgrid` everywhere.** When plotting with `contourf` or
`pcolor` — which expect meshgrid orientation — transpose: `contourf(xc, yc, PHI.')`.

---

## 3.6 Reading the sparsity pattern

```matlab
S = poisson2d(8, 6, 1/8, 1/6, {'D','D','N','N'});
figure; spy(S.A); title('2D Laplacian sparsity');
fprintf('size %d x %d, %d nonzeros (%.2f%% dense)\n', ...
        size(S.A,1), size(S.A,2), nnz(S.A), 100*nnz(S.A)/numel(S.A));
```

You should see five diagonals: the main one, $\pm 1$ (x-neighbours), and
$\pm n_x$ (y-neighbours). Gaps in the $\pm 1$ diagonals occur exactly at the
x-boundaries — cell $(n_x, j)$ must not connect to cell $(1, j+1)$ even though
they are adjacent in the flattened vector. **`kron` handles this correctly and
hand-rolled loops usually don't.** If you ever write the assembly by hand, this
"wraparound" is the bug to look for; its signature is a solution with a spurious
diagonal streak.

---

## Exercises

**3.1** Run `poisson_mms`. Confirm second order. Then change `ndgrid` to
`meshgrid` and rerun. Explain what you see. Then set `ny = nx` and rerun the
`meshgrid` version — does the bug still show?

**3.2** Solve steady conduction in a 200 µm × 100 µm chip cross-section: hot wall
$T=350$ K at the bottom, $T=300$ K at the top, insulated ($\partial T/\partial n = 0$)
sides, no volumetric source. Verify the solution is a pure linear gradient in
$y$, independent of $x$. Compute the heat flux through the top wall and compare
to $k\Delta T/H$.

**3.3** Add a uniform heat source $q$ (a Joule-heated resistor) to 3.2. The
analytic solution is now a parabola. Verify against it.

**3.4 (the singular one)** Change all four BCs in 3.2 to Neumann. Comment out the
pinning in `poisson2d` and try to solve. Record the exact warning MATLAB gives.
Restore the pinning, solve again, and confirm that although $\phi$ changes when
you pin a different cell, $\nabla\phi$ does not. *(This is the fact that makes
pinning legitimate — verify it, don't take my word for it.)*

**3.5** Time the assembly vs the solve:
```matlab
n=128; tic, S=poisson2d(n,n,1/n,1/n,{'D','D','D','D'}); t_asm=toc;
f=randn(n,n); z={0,0,0,0};
tic, for k=1:100, phi=S.solve(f,z); end, t_slv=toc/100;
fprintf('assembly %.3f s, solve %.4f s, ratio %.0f\n',t_asm,t_slv,t_asm/t_slv);
```
How many timesteps before factorizing-once has paid for itself?

**3.6** Verify the compatibility condition. Set up an all-Neumann problem with
$f = 1$ everywhere and $g = 0$ on all walls. Note $\int f \ne \oint g$. Solve it
anyway and look at the result. Then set $f = 1 - \langle f\rangle$ (mean-zero)
and solve again. *(This is exactly the failure mode of a badly-posed outflow BC
in Tutorial 8.)*

---

**Next:** [Tutorial 4 — Unsteady diffusion and stability](04-diffusion-and-stability.md).
You can solve steady problems. Now add time — and discover that the timestep is
not yours to choose freely.
