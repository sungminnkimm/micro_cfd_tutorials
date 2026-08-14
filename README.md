# Microscale CFD in MATLAB — A Build-It-Yourself Course

You know fluid mechanics. You have never written a solver. By the end of this
course you will have written, from a blank editor, an incompressible
Navier–Stokes solver on a staggered grid, coupled it to scalar transport and to
a level-set interface tracker, and used it to answer real microfluidic design
questions — droplet size in a T-junction, mixing length in a serpentine channel,
dispersion in a capillary.

Everything runs in base MATLAB. No toolboxes required.

---

## How to use this course

Each tutorial is a markdown file in `tutorials/`. Each has matching runnable
code in `code/chXX/`.

**Read the tutorial with MATLAB open.** Type the code. Do not copy-paste the
short snippets — typing them is how the indexing conventions get into your
fingers, and indexing conventions are 80% of what makes CFD code hard. The long
solvers you may copy; they are provided complete in `code/`.

Every tutorial ends with **Exercises**. They are not optional garnish. The
progression from "I ran the demo" to "I can build my own model" happens entirely
inside those exercises.

### Setup

```matlab
cd D:\Projects\CFD\code
setup_path            % adds lib/ and all chapter folders to the MATLAB path
```

Run it once per MATLAB session. Then any script can be launched by name:

```matlab
cavity
poiseuille
mms_driver
```

### Regression test

```matlab
run_all_tests                 % all 30 scripts
MODE='fast'; run_all_tests    % the 22 quick ones (~3 min)
```

Run this after modifying anything in `lib/`. The quick set covers every
verification check; the slow set adds the benchmark and two-phase demos.

---

## Verified results

Every number below is actual output from the code in this repository, not an
aspiration. Reproduce them with `run_all_tests`.

| Test | Expected | Measured |
|---|---|---|
| FD stencils (T2) | 1, 2, 2 | **1.00, 2.00, 2.00** |
| Poisson MMS, non-square grid (T3) | 2 | **1.99** |
| Time integration FE/BE/CN (T4) | 1, 1, 2 | **1.00, 1.00, 2.00** |
| Upwind numerical diffusivity (T5) | $u\Delta x(1{-}C)/2$ | **6.250e-04 vs 6.250e-04** |
| Poiseuille order (T6) | 2 | **2.00**, $u_{max}/\bar u$ = 1.498 |
| Cubic law $Q\propto h^3$ (T6) | 3 | **3.0000** |
| Cavity vs Ghia $Re$=100 (T7) | — | **err < 0.005**, `divMax` 7e-15 |
| **NS solver MMS (T13)** | **2** | **2.00 (u and v)**, `divMax` 1e-15 |
| Channel $\Delta p$ (T8) | $12\mu\bar uL/h^2$ | **within 0.09%**, `divMax` 4e-14 |
| Taylor–Aris $D_{eff}$ (T9) | $D(1+2Pe^2/105)$ | **within 1.4%** |
| Static droplet $\Delta p$ (T11) | $\sigma/R$ | **within 1.0%** at 64² |

The MMS row is the strongest: a first-order bug anywhere in the solver would
drag that 2.00 to 1.00.

### And two results that do *not* pass

Stated up front, because they determine how Part III should be read.

**1. Spurious currents.** `static_droplet.m` measures $5.8\times10^{-3}$ m/s at
64² for the default water-in-oil properties — **six times larger** than the 1
mm/s flow a microchannel simulation would try to resolve. A real limitation of
the level-set/CSF implementation; Tutorial 11.4 gives the options.

**2. The T-junction produces no detached droplets** at the affordable
resolution. The run is numerically healthy (`divMax` ≈ 2e-12, dispersed area
growing as injection demands), but at 16.7 cells across the channel the neck is
too thin to resolve pinch-off, and §12.4 asks for ≥20 *across the neck*.

The second one is deliberate rather than papered over. Had it pinched off at this
grid, breakup would have been triggered by two interfaces coming within a cell of
each other — a **numerical** criterion, not a physical one — and the droplet size
would have been a function of $\Delta x$. **A demo that silently produced a
plausible droplet length here would be worse than one that refuses to.**
Exercise 12.5 is where you refine until it works and the answer stops moving.

**So: trust Part III for regimes, trends, and mechanism. Not for absolute
droplet volumes.**

### A note on the `QUIET` flag

Every script begins with

```matlab
if ~exist('QUIET','var'), QUIET = false; end
```

Plotting and animation are wrapped in `if ~QUIET`. Set `QUIET = true` before
running if you want numbers only (useful in loops, batch runs, and convergence
studies). Leave it alone for normal interactive use.

---

## Curriculum

### Part I — Learning to discretize (Tutorials 1–5)

You cannot debug a Navier–Stokes solver if you cannot debug a Poisson solver.
Part I builds the numerical machinery on problems simple enough that you always
know the right answer.

| # | Tutorial | What you build | Key idea |
|---|----------|----------------|----------|
| 1 | [Why microscale is different](tutorials/01-why-microscale-is-different.md) | A scaling calculator | Re, Pe, Ca, Bo — which terms you are allowed to drop |
| 2 | [Finite differences and truncation error](tutorials/02-finite-differences.md) | Derivative stencils + a convergence study | If your error doesn't drop at the right rate, your code is wrong |
| 3 | [Boundary conditions and sparse linear systems](tutorials/03-poisson-and-sparse-systems.md) | A 2D Poisson solver | Assembling `A\b`; Dirichlet vs Neumann; the singular-matrix trap |
| 4 | [Unsteady diffusion and stability](tutorials/04-diffusion-and-stability.md) | Explicit, implicit, Crank–Nicolson | Von Neumann analysis; why your timestep blew up |
| 5 | [Advection and the Péclet number](tutorials/05-advection-and-peclet.md) | Advection–diffusion solver | Upwinding, numerical diffusion, cell Péclet, the wiggles |

### Part II — The Navier–Stokes solver (Tutorials 6–8)

| # | Tutorial | What you build | Key idea |
|---|----------|----------------|----------|
| 6 | [Stokes flow and the staggered grid](tutorials/06-stokes-and-staggered-grid.md) | Poiseuille solver, analytic benchmarks | Why colocated grids checkerboard; the MAC grid |
| 7 | [The projection method](tutorials/07-projection-method.md) | **Full incompressible NS solver** | Fractional step; the pressure Poisson equation; lid-driven cavity validation |
| 8 | [Building your own geometry](tutorials/08-microchannel-geometry.md) | Channel, contraction, obstacle, serpentine | Inlet/outlet BCs; masking solid regions; validating against Poiseuille |

### Part III — Microscale phenomena (Tutorials 9–12)

| # | Tutorial | What you build | Key idea |
|---|----------|----------------|----------|
| 9 | [Mixing and dispersion](tutorials/09-mixing-and-dispersion.md) | Passive scalar on your flow field | Why microchannels don't mix; Taylor–Aris; mixing index; chaotic advection |
| 10 | [Surface tension and capillarity](tutorials/10-surface-tension-capillarity.md) | Young–Laplace solver, meniscus shapes | Curvature, contact angle, capillary pressure, Ca and Bo |
| 11 | [Two-phase flow with level sets](tutorials/11-level-set-two-phase.md) | **Interface tracker + CSF surface tension** | Signed distance, reinitialization, spurious currents, the capillary timestep |
| 12 | [Droplet generation](tutorials/12-droplet-generation.md) | T-junction and flow-focusing | Squeezing / dripping / jetting; droplet-size scaling laws |

### Part IV — Doing it yourself (Tutorials 13–14)

| # | Tutorial | What you build | Key idea |
|---|----------|----------------|----------|
| 13 | [Verification and validation](tutorials/13-verification-and-validation.md) | MMS driver, grid convergence index | How to prove your own code is right |
| 14 | [From a question to a model](tutorials/14-your-own-model.md) | Your project | Non-dimensionalization, resolution budgets, a debugging playbook |

Tutorial 13 closes with **three real bugs caught while building this course** —
a missing mask in the pressure operator (`divMax` = 416), an uncorrected outflow
face (`divMax` = 1.5), and a flipped curvature sign (pressure jump of the wrong
sign). All three produced plausible-looking output; two gave a *correct* answer
to the first quantity anyone would check. They were caught only by diagnostics
that print a number with a known right answer, which is the habit this course is
really teaching.

### Appendices

- [A. MATLAB performance for CFD](tutorials/A-matlab-performance.md) — vectorization, sparse patterns, preallocation, profiling
- [B. Common bugs and their symptoms](tutorials/B-common-bugs.md) — a lookup table from "what I see" to "what's wrong"
- [C. Notation and symbols](tutorials/C-notation.md)

---

## The physical setting

Unless stated otherwise, "microscale" here means:

| Quantity | Symbol | Typical value |
|----------|--------|---------------|
| Channel width/depth | $L$ | 10–500 µm |
| Mean velocity | $U$ | 0.1–100 mm/s |
| Fluid | water | $\rho = 998$ kg/m³, $\mu = 1.0\times10^{-3}$ Pa·s |
| Small-molecule diffusivity | $D$ | $\sim 5\times10^{-10}$ m²/s |
| Water/air surface tension | $\sigma$ | 0.072 N/m |

which puts you at $Re \sim 10^{-3}$–$10$, $Pe \sim 10$–$10^4$, $Ca \sim 10^{-5}$–$10^{-2}$.
Tutorial 1 explains what each of those numbers buys you.

---

## What this course does not cover

Deliberate omissions, so you know where the edges are:

- **Electrokinetics** (electroosmotic flow, EDL, dielectrophoresis). Big topic,
  genuinely important in microfluidics, and orthogonal to what's here. Tutorial 14
  sketches how to bolt a Poisson–Boltzmann solve onto the Tutorial 7 solver.
- **Rarefied / slip flow** (Knudsen number, gas microflows). Only matters for gases
  below ~1 µm. Same note in Tutorial 14.
- **3D.** Everything here is 2D (and 2D-axisymmetric where noted). The numerics
  are identical in 3D; the cost is not. Tutorial 14 discusses when 2D lies to you
  — and for microchannels with low aspect ratio, it lies quite a lot.
- **Unstructured meshes and finite volume on arbitrary polygons.** We use
  structured Cartesian grids with finite differences / finite volumes. This is the
  right choice for learning and for most microchannel geometries.
- **Turbulence.** At $Re \sim 1$ there isn't any. This is the one thing microscale
  CFD makes *easier*.

