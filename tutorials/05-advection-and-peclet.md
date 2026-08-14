# Tutorial 5 — Advection and the Péclet Number

**Goal:** Transport a scalar with the flow. Discover why the obvious scheme
fails, why the obvious fix is also bad, and how to tell numerical diffusion from
the real thing. At $Pe = 200$ this is the tutorial that decides whether your
mixing predictions mean anything.

**Code:** `code/ch05/advect1d.m`, `code/ch05/cell_peclet.m`, `code/ch05/num_diffusion.m`

---

## 5.1 The equation and why it's harder

$$\frac{\partial c}{\partial t} + u\frac{\partial c}{\partial x} = D\frac{\partial^2 c}{\partial x^2}$$

Diffusion (Tutorial 4) is **parabolic**: it smooths, it forgives, errors decay.
Advection is **hyperbolic**: it translates without smoothing, and errors *persist
and propagate*. There is no dissipation to hide your mistakes.

The pure-advection exact solution is trivial — $c(x,t) = c_0(x - ut)$, the initial
profile slid along, unchanged. That "unchanged" is exactly what your scheme will
fail to reproduce.

---

## 5.2 Attempt 1: central differencing (and why it fails)

Reuse the second-order central stencil from Tutorial 2:

$$c_i^{n+1} = c_i^n - \frac{u\Delta t}{2\Delta x}\left(c_{i+1}^n - c_{i-1}^n\right)$$

Von Neumann analysis, with the Courant number $C = u\Delta t/\Delta x$:

$$G = 1 - \frac{C}{2}\left(e^{ik\Delta x} - e^{-ik\Delta x}\right) = 1 - iC\sin(k\Delta x)$$

$$|G|^2 = 1 + C^2\sin^2(k\Delta x) \ \ge\ 1$$

**$|G| > 1$ for every wavenumber and every timestep.** Unconditionally *unstable*.
The second-order-accurate scheme is useless on its own. This is the first time
in the course that "more accurate" and "works" come apart, and it's worth
sitting with: accuracy and stability are independent properties.

(With enough physical diffusion present, the diffusion term's damping can
stabilize it — that's §5.4. But central differencing for pure advection is dead.)

---

## 5.3 Attempt 2: upwind differencing

The physical insight: information in advection travels *one way*, downstream.
So take the derivative from the upwind side:

$$\frac{\partial c}{\partial x}\bigg|_i \approx
\begin{cases}
\dfrac{c_i - c_{i-1}}{\Delta x} & u > 0 \quad\text{(look left)}\\[2ex]
\dfrac{c_{i+1} - c_i}{\Delta x} & u < 0 \quad\text{(look right)}
\end{cases}$$

Vectorized in MATLAB without branching:

```matlab
up = max(u,0);  um = min(u,0);
dcdx = up.*(c - circshift(c,1))/dx + um.*(circshift(c,-1) - c)/dx;
```

Von Neumann for $u>0$:

$$G = 1 - C\left(1 - e^{-ik\Delta x}\right)$$
$$|G|^2 = 1 - 2C(1-C)(1-\cos k\Delta x)$$

$|G|\le 1$ requires $C(1-C) \ge 0$, i.e.

$$\boxed{C = \frac{u\Delta t}{\Delta x} \le 1} \qquad \text{(the CFL condition)}$$

**Physical meaning: a fluid parcel may not cross more than one cell per
timestep.** If it does, the scheme is trying to get information from a cell it
never looked at. Courant–Friedrichs–Lewy, 1928, and it constrains every explicit
advection scheme ever written.

### The catch: upwinding is diffusion in disguise

Upwinding is stable. It is also **first-order** accurate, and the error it makes
is not random — it is specifically diffusive. Taylor-expand the upwind stencil:

$$\frac{c_i - c_{i-1}}{\Delta x} = \frac{\partial c}{\partial x} - \frac{\Delta x}{2}\frac{\partial^2 c}{\partial x^2} + O(\Delta x^2)$$

Substituting back, the scheme you are *actually* solving is

$$\frac{\partial c}{\partial t} + u\frac{\partial c}{\partial x} = \left(D + \underbrace{\frac{u\Delta x}{2}(1-C)}_{D_{num}}\right)\frac{\partial^2 c}{\partial x^2}$$

You are solving the right equation with the **wrong diffusivity**. The extra
piece $D_{num}$ is **numerical diffusion**, and it is the single most dangerous
error in microfluidic CFD.

### Why it's dangerous: run the numbers

Microchannel: $u = 1$ mm/s, $\Delta x = 1$ µm, $D = 5\times10^{-10}$ m²/s.

$$D_{num} \approx \frac{u\Delta x}{2} = \frac{10^{-3}\times10^{-6}}{2} = 5\times10^{-10}\ \text{m}^2\text{/s}$$

$$\frac{D_{num}}{D} \approx 1$$

**Your numerical diffusion equals the real thing.** Your simulation will show a
device mixing twice as fast as it does. And the failure is silent: the result
looks smooth, plausible, and publishable.

A useful way to say the same thing: the ratio $D_{num}/D$ is just the **cell
Péclet number over 2**:

$$\frac{D_{num}}{D} = \frac{u\Delta x}{2D} = \frac{Pe_{cell}}{2}$$

So the requirement "numerical diffusion below 10% of physical" means
$Pe_{cell} < 0.2$, which at $Pe = 200$ over a 100 µm channel demands
$\Delta x < 0.1$ µm — **1000 cells across a single channel.** That is why
first-order upwinding is not good enough for microfluidic mixing, and why §5.5
exists.

---

## 5.4 The cell Péclet number and the wiggle problem

For the steady advection–diffusion equation with central differencing, the
discrete solution is oscillation-free only if

$$Pe_{cell} = \frac{|u|\Delta x}{D} \le 2$$

Above that, central differencing produces spurious cell-to-cell oscillations —
**wiggles** — near sharp gradients. They are not instability (the solution
doesn't blow up); they are the scheme's inability to represent a boundary layer
thinner than a cell.

This puts you in a genuine bind:

| Scheme | Wiggles? | Numerical diffusion? |
|---|---|---|
| Central | Yes, if $Pe_{cell} > 2$ | None |
| Upwind | Never | $D_{num} = u\Delta x/2$ — huge at high $Pe$ |

**Neither is acceptable at $Pe = 200$.** The way out is a scheme that is
high-order where the solution is smooth and drops to upwind only near sharp
gradients. That's §5.5.

`code/ch05/cell_peclet.m` demonstrates the bind on the classic 1D boundary-layer
problem, which has the exact solution

$$c(x) = \frac{1 - e^{Pe\,x/L}}{1 - e^{Pe}}$$

— a flat region with an exponential boundary layer of thickness $L/Pe$ at the
outflow end. Run it at $Pe_{cell} = 1, 2, 4, 20$ and watch central differencing
disintegrate while upwind stays smooth but smears the layer.

---

## 5.5 Doing it properly: flux limiters (TVD)

The fix is to write the scheme in **flux form** and blend two fluxes based on how
smooth the solution locally is.

Conservative flux form — mass in = mass out, exactly:

$$c_i^{n+1} = c_i^n - \frac{\Delta t}{\Delta x}\left(F_{i+1/2} - F_{i-1/2}\right)$$

Any scheme written this way conserves total mass to machine precision, which is
worth having for its own sake.

Now build the face flux as low-order plus a limited correction:

$$F_{i+1/2} = \underbrace{u\,c_i}_{\text{upwind}} + \underbrace{\tfrac{1}{2}u(1-C)\,\psi(r_i)\,(c_{i+1}-c_i)}_{\text{anti-diffusive correction}}$$

where $r_i$ is the **ratio of consecutive gradients**,

$$r_i = \frac{c_i - c_{i-1}}{c_{i+1} - c_i}$$

which measures local smoothness: $r\approx 1$ means smooth, $r\le 0$ means an
extremum. The **limiter** $\psi(r)$ decides how much correction to allow:

| Limiter | $\psi(r)$ | Character |
|---|---|---|
| Upwind | $0$ | most diffusive |
| Lax–Wendroff | $1$ | 2nd order, oscillatory |
| **van Leer** | $\dfrac{r+\|r\|}{1+\|r\|}$ | smooth, good default |
| **Superbee** | $\max(0,\min(2r,1),\min(r,2))$ | sharpest, can over-steepen |
| minmod | $\max(0,\min(r,1))$ | most cautious |

All satisfy $\psi(r) = 0$ for $r\le0$ (kill the correction at extrema — no new
maxima can be created, which is the **TVD** property, hence no wiggles) and
$\psi(1) = 1$ (full second-order accuracy where the solution is smooth).

```matlab
function psi = limiter(r, name)
switch lower(name)
    case 'upwind',   psi = zeros(size(r));
    case 'lw',       psi = ones(size(r));
    case 'vanleer',  psi = (r + abs(r))./(1 + abs(r));
    case 'superbee', psi = max(0, max(min(2*r,1), min(r,2)));
    case 'minmod',   psi = max(0, min(r,1));
end
psi(~isfinite(psi)) = 0;      % r = Inf/NaN where the denominator vanishes
end
```

That last line matters. When $c_{i+1} = c_i$ the ratio is $0/0$ or $x/0$. Guard
it — an unguarded `NaN` propagates through the whole field in one step and is a
classic first-run failure.

**Result:** second-order accurate in smooth regions, no oscillations at
discontinuities, and *dramatically* less numerical diffusion than upwind.
`advect1d.m` measures this: for a top-hat advected once around a periodic
domain, van Leer typically preserves several times the peak amplitude that
first-order upwind does, at identical cost.

---

## 5.6 Measuring your own numerical diffusion

You should never *assume* your numerical diffusion is acceptable. Measure it.

**The method:** advect a Gaussian of known initial width $\sigma_0$ with $D = 0$
(no physical diffusion). Exactly, it should translate with $\sigma$ unchanged.
Numerically it spreads. Since a diffusing Gaussian satisfies
$\sigma^2 = \sigma_0^2 + 2Dt$, you can read off the effective diffusivity:

$$D_{num} = \frac{\sigma^2(t) - \sigma_0^2}{2t}$$

`code/ch05/num_diffusion.m` does exactly this, comparing the measured $D_{num}$
against the theoretical $u\Delta x(1-C)/2$ for upwind, and reporting the ratio
$D_{num}/D_{phys}$ for a realistic microchannel.

**Do this for every mixing simulation you ever run.** If $D_{num}/D_{phys} > 0.1$,
your mixing result is dominated by your grid, not your physics. This one
diagnostic separates trustworthy microfluidic CFD from the other kind.

### Reporting it honestly

The right sentence for a report reads: *"With $\Delta x = 0.5$ µm and van Leer
limiting, the measured numerical diffusivity was $3\times10^{-11}$ m²/s, 6% of the
physical value $5\times10^{-10}$ m²/s; halving the grid changed the predicted
mixing length by 2%."* That is a defensible claim. "We used a fine grid" is not.

---

## 5.7 Timestep: now you have two limits

Advection–diffusion explicit stability needs **both**:

$$\Delta t \le \frac{\Delta x}{|u|_{max}} \quad\text{(CFL)} \qquad\text{and}\qquad \Delta t \le \frac{\Delta x^2}{2D}\ \text{(1D)}$$

Which binds depends on $Pe_{cell}$:

- $Pe_{cell} < 2$: diffusion limit binds. Common for **momentum** at low $Re$.
- $Pe_{cell} > 2$: CFL binds. Common for **scalars** at high $Pe$.

Take the minimum, with a safety factor:

```matlab
dt = 0.4 * min( dx/max(abs(u(:))), dx^2/(2*D) );
```

The 0.4 is not superstition — it covers the fact that the two limits interact,
and that in 2D each is tighter than its 1D form. Appendix B lists the full set of
timestep limits used in this course, including the capillary one from Tutorial 11.

---

## Exercises

**5.1** Run `advect1d` with `'central'`. Watch it blow up. Confirm the growth
rate matches $|G| = \sqrt{1 + C^2\sin^2(k\Delta x)}$ for the sawtooth mode.

**5.2** Run `advect1d` with upwind at $C = 0.5$ and $C = 0.95$. Which preserves
the top hat better? Explain using the $(1-C)$ factor in $D_{num}$. *(The
counter-intuitive answer — bigger timestep, less numerical diffusion — is real,
and at $C=1$ upwind is exact. Understand why before relying on it.)*

**5.3** Run all limiters on the same problem. Rank them by peak preservation and
by overshoot. Confirm superbee is sharpest and that it slightly squares off a
smooth Gaussian.

**5.4** Run `cell_peclet` at $Pe_{cell} = 1, 2, 4, 20$. Find the exact threshold
where central differencing starts to wiggle.

**5.5 (the one that matters)** Take your Tutorial 1 device. Compute
$Pe_{cell} = u\Delta x/D$ for a grid of 50 cells across the channel. Then compute
$D_{num}/D$ for upwind. Now find the $\Delta x$ that gets $D_{num}/D$ below 0.1
with upwind, and again with van Leer. How many cells across the channel does each
need? *(This calculation is what should size every mixing grid you ever build.)*

**5.6** Modify `num_diffusion` to measure $D_{num}$ for the van Leer scheme.
It is not a constant — it depends on the local solution smoothness, so a Gaussian
and a top hat give different answers. What does that imply about quoting a single
number for a nonlinear scheme?

---

**Next:** [Tutorial 6 — Stokes flow and the staggered grid](06-stokes-and-staggered-grid.md).
You now have every piece needed for the momentum equations. Part II assembles
them into a real Navier–Stokes solver.
