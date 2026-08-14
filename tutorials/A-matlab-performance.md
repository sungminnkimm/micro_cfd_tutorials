# Appendix A — MATLAB Performance for CFD

MATLAB is fast enough for everything in this course *if* you avoid a handful of
specific mistakes. Here they are, in order of impact.

---

## A.1 Factorize once (the big one)

Measured in this course: the pressure Poisson solve is called once per timestep
with the **same matrix** and a different right-hand side. Refactorizing every
step is the single most expensive mistake you can make.

```matlab
% BAD  -- refactorizes every call, O(n^1.5) each time
for n = 1:nsteps
    phi = A \ b;
end

% GOOD -- factorize once
dA = decomposition(A, 'lu');       % or [L,U,P,Q,R] = lu(A)
for n = 1:nsteps
    phi = dA \ b;
end
```

`decomposition` (R2017b+) is the clean modern form and it caches the
factorization. For a 128² Poisson operator this is a 10–50× speedup on the whole
simulation.

**This only works if the matrix is constant.** Adaptive timestepping with an
implicit viscous term forces refactorization every time $\Delta t$ changes,
which usually costs more than the adaptivity saves. Pick a fixed $\Delta t$ —
this is why `ns_solver` does.

---

## A.2 Vectorize the stencils

```matlab
% BAD -- ~100x slower
for i = 2:nx-1
    for j = 2:ny-1
        lap(i,j) = (c(i+1,j)-2*c(i,j)+c(i-1,j))/dx^2 + ...
    end
end

% GOOD
lap = (c(3:end,2:end-1) - 2*c(2:end-1,2:end-1) + c(1:end-2,2:end-1))/dx^2 ...
    + (c(2:end-1,3:end) - 2*c(2:end-1,2:end-1) + c(2:end-1,1:end-2))/dy^2;
```

Array slicing is the idiom the whole course uses. It is worth taking the time to
read a slice expression fluently — `c(3:end,2:end-1)` is "shift +1 in x, interior
in y."

MATLAB's JIT has narrowed the loop/vectorized gap considerably in recent
versions, but slicing is also *clearer* for stencils once you're used to it, and
it maps directly onto the mathematical expression.

---

## A.3 Preallocate

```matlab
% BAD -- reallocates and copies every iteration
hist = [];
for n = 1:nsteps, hist(end+1) = ke; end

% GOOD
hist = zeros(nsteps,1);
for n = 1:nsteps, hist(n) = ke; end
```

MATLAB flags this in the editor. Heed it — for a 100k-step droplet run, growing
an array is genuinely significant.

---

## A.4 Sparse matrices

```matlab
A = zeros(n);        % 128^2 grid -> 16384^2 doubles = 2 GB. Don't.
A = sparse(n,n);     % stores only nonzeros
```

Build with **triplets**, not element-by-element assignment:

```matlab
% BAD -- each assignment reallocates the sparse structure
A = sparse(n,n);
for k = 1:n, A(k,k) = -4; end

% GOOD -- one allocation
I = zeros(5*n,1); J = I; V = I; ptr = 0;
% ... fill I, J, V ...
A = sparse(I(1:ptr), J(1:ptr), V(1:ptr), n, n);
```

`poisson2d.m` uses the triplet form for exactly this reason (and because the
masked variant can't be expressed as a clean `kron`).

---

## A.5 Profiling: measure, don't guess

```matlab
profile on
cavity
profile viewer
```

Or time sections directly:

```matlab
t = zeros(1,3);
for n = 1:nsteps
    a=tic; predictor;  t(1)=t(1)+toc(a);
    a=tic; poisson;    t(2)=t(2)+toc(a);
    a=tic; corrector;  t(3)=t(3)+toc(a);
end
fprintf('predictor %.1f%%  poisson %.1f%%  corrector %.1f%%\n', 100*t/sum(t));
```

**Typical split for this course's solver:** the Poisson solve dominates at
60–80%, which is why A.1 matters so much and why everything else is secondary.

---

## A.6 When direct solves run out

Direct LU is fine up to ~256². Beyond that, memory for the factors becomes the
problem (fill-in grows faster than the matrix).

**Iterative with a preconditioner:**

```matlab
Lic = ichol(-A);                       % incomplete Cholesky (A is negative definite)
[phi, flag] = pcg(-A, -b, 1e-10, 200, Lic, Lic', phiPrev);
```

Passing `phiPrev` as the initial guess matters enormously in a time loop — the
pressure changes little between steps, so you typically converge in a handful of
iterations instead of a hundred.

**FFT/DCT-based fast Poisson solver.** For a uniform grid with uniform BCs, the
Laplacian is diagonalized by the discrete cosine transform:

```matlab
% Homogeneous Neumann on all sides
lam = -4*(sin(pi*(0:nx-1)'/(2*nx)).^2/dx^2 + sin(pi*(0:ny-1)/(2*ny)).^2/dy^2);
lam(1,1) = 1;                          % the null mode
phi = idct2( dct2(b) ./ lam );
```

$O(N\log N)$ and very fast — but it only works for constant coefficients and
simple BCs, so it can't handle the solid mask from Tutorial 8.

---

## A.7 Memory

For an $n_x\times n_y$ grid, a double array is $8n_xn_y$ bytes. A 512² field is
2 MB — trivial. Storing 1000 timesteps of it is 2 GB — not trivial.

**Save selectively:**

```matlab
saveEvery = 100;
nSave = floor(nsteps/saveEvery);
snap = zeros(nx, ny, nSave);
for n = 1:nsteps
    ...
    if mod(n,saveEvery)==0, snap(:,:,n/saveEvery) = c; end
end
```

Or use `single` for output arrays you'll only plot — halves the memory and
plotting can't tell.

---

## A.8 Things that do *not* help

Worth stating, since they absorb time that A.1 deserves:

- **`parfor` on the time loop.** Timesteps are sequential; you cannot parallelize
  them. (Parameter *sweeps* are a different matter and parallelize perfectly.)
- **GPU arrays** for problems this size. Transfer overhead dominates below
  ~1000². `gpuArray` pays off for large 3D runs, not 2D microchannels.
- **Micro-optimizing the predictor** when the Poisson solve is 70% of runtime.
  Profile first.
- **Converting to C/MEX** before you have verified the MATLAB version. Get it
  right, then get it fast — and keep the slow version as a reference to test
  against.

---

## A.9 A realistic performance budget

Measured on the machine used to build this course:

| Case | Grid | Steps | Time |
|---|---|---|---|
| Lid-driven cavity, $Re=100$ | 64² | 5,400 | 3 s |
| Lid-driven cavity, $Re=100$ | 128² | 21,600 | 45 s |
| Straight channel | 160×40 | 1,900 | 2 s |
| MMS verification (5 grids) | to 64² | 132,000 | ~90 s |
| Static droplet (two-phase) | 64² | 30,000 | 41 s |

Note the cavity scaling: 4× the cells and 4× the steps gives ~16× the time — the
$N^4$ behaviour Tutorial 4 predicted for explicit diffusion. That table is the
one to consult when deciding whether a planned run is affordable, and
`regime_map.m` does the same arithmetic for droplet problems before you commit.
