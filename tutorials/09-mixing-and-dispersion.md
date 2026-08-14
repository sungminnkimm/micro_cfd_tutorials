# Tutorial 9 — Mixing and Dispersion

**Goal:** Answer the question microfluidics actually asks: *how long must my
channel be for these two streams to mix?* Learn why the answer is
disappointing, and what geometries fix it.

**Code:** `code/ch09/mixing_study.m`, `code/ch09/taylor_dispersion.m`, `code/lib/scalar_step.m`

---

## 9.1 Why microchannels refuse to mix

Three facts from Part I, combined:

1. $Re \ll 1$ → **no turbulence.** The only mixing mechanism available to you
   at the macro scale is gone.
2. $Pe \gg 1$ → **advection dominates diffusion.** Streams flow side by side.
3. Stokes flow is **reversible** → stirring can be undone. Simply moving fluid
   around does not, by itself, mix.

The consequence, from Tutorial 1: two streams in a 100 µm channel at 1 mm/s need
**~20 mm** of channel to mix by diffusion alone. Tutorial 4's more careful
calculation refined that to $\approx 0.3\,Pe\,w$, but it's the same story — a
significant fraction of your entire chip, spent on mixing.

And it gets worse for large molecules. $D \propto 1/R$ (Stokes–Einstein), so:

| Species | $D$ [m²/s] | $Pe$ | $L_{mix}$ (100 µm, 1 mm/s) |
|---|---|---|---|
| Small ion | $2\times10^{-9}$ | 50 | 1.5 mm |
| Small molecule | $5\times10^{-10}$ | 200 | 6 mm |
| Protein | $4\times10^{-11}$ | 2500 | 75 mm |
| 1 µm bead | $4\times10^{-13}$ | $2.5\times10^5$ | **7.5 m** |

That last row is the design constraint that created the entire field of active
and chaotic micromixers.

---

## 9.2 The equation, and the shortcut that makes it affordable

$$\frac{\partial c}{\partial t} + \mathbf{u}\cdot\nabla c = D\nabla^2 c$$

with $\mathbf{u}$ from your Tutorial 7/8 solver.

**The shortcut:** at $Re \ll 1$ with steady boundary conditions, the velocity
field is *constant in time*. So:

1. Solve Navier–Stokes **once** to steady state. Save `u`, `v`.
2. Transport the scalar on that frozen field, with a timestep chosen for the
   *scalar* equation only.

This is not an approximation — it's exact for steady Stokes flow. And it is a
huge saving, because as Tutorial 7 noted, the viscous timestep limit is $Sc\approx2000$
times more restrictive than the scalar one. Freezing the flow buys you back that
factor.

`code/ch08/serpentine.m` saves its flow field to `serpentine_flow.mat` for
exactly this reason.

---

## 9.3 Numerics: this is where results are won or lost

From Tutorial 5: the scalar equation at $Pe = 200$ is advection-dominated, and
**first-order upwinding produces numerical diffusion comparable to the physical
diffusivity.** A mixer simulated with upwind on a coarse grid will appear to work
beautifully and be entirely wrong.

The measured numbers from Tutorial 5's `num_diffusion.m`:

| Scheme | Measured $D_{num}$ (that test) | Relative |
|---|---|---|
| upwind | $6.25\times10^{-4}$ | 1× |
| minmod | $1.83\times10^{-5}$ | 1/34 |
| **van Leer** | $8.26\times10^{-7}$ | **1/757** |
| superbee | $-1.69\times10^{-5}$ | anti-diffusive |

**Use van Leer.** Same cost per step as upwind, ~750× less numerical diffusion.

And **always report your measured $D_{num}/D_{phys}$.** `mixing_study.m`
computes it for the grid you actually used. If it exceeds ~0.1, your mixing
result is a property of your mesh, not your device.

### Timestep for the scalar

$$\Delta t \le 0.4\min\left(\frac{\Delta x}{|u|_{max}},\ \frac{1}{2D}\left(\frac{1}{\Delta x^2}+\frac{1}{\Delta y^2}\right)^{-1}\right)$$

At high $Pe$ the CFL term binds — which is convenient, because it's the gentler
of the two.

---

## 9.4 Measuring mixing

### The mixing index

$$M = \frac{\sigma(c)}{\sigma(c_0)}, \qquad \sigma = \text{standard deviation}$$

$M=1$ unmixed, $M=0$ perfect. Usually evaluated cross-section by cross-section
down the channel.

**A caveat that matters.** $M$ is blind to *structure*. A field striped into
many thin lamellae of pure A and pure B has the same $\sigma$ as one thick
stripe of each — but it is about to mix vastly faster, because diffusion time
scales as (stripe width)². So when a mixer works by *stretching and folding*,
$M$ will understate its effect until diffusion finally cashes in the extra
interface.

Therefore also track the **interfacial length** $\int|\nabla c|\,dA$. Together,
$M$ and interfacial length tell you both *how mixed* you are and *how fast you
are about to become mixed*. `mixing_study.m` reports both.

### Beware the convention

Some papers report **mixing efficiency** $1-M$. So "95% mixed" can mean $M=0.05$
or $M=0.95$ depending on the author. State which you mean. This course always
reports $M$ as defined above.

---

## 9.5 Taylor–Aris dispersion: the other thing that spreads your sample

Separate phenomenon from cross-stream mixing, and just as important.

A plug of solute injected into Poiseuille flow spreads **along** the channel far
faster than molecular diffusion alone would explain. Why: fluid at the
centreline moves at $1.5\bar u$ while fluid near the walls barely moves.
Molecules sample different streamlines by diffusing across the channel, and the
combination of shear + cross-stream diffusion produces an *effective* axial
diffusivity:

$$\boxed{D_{eff} = D\left(1 + \frac{2}{105}Pe_h^2\right)}, \qquad Pe_h = \frac{\bar u h}{D}$$

for flow between parallel plates separated by $2h$. (For a circular tube the
coefficient is $1/48$ instead of $2/105$.)

**Run the numbers.** $\bar u = 1$ mm/s, $h = 50$ µm, $D = 5\times10^{-10}$:

$$Pe_h = 100, \qquad D_{eff} = D(1 + 190) = 191 D$$

Your sample plug spreads **191 times faster** than molecular diffusion. For any
device that injects a discrete plug — chromatography, electrophoresis, droplet
trains, sequential sample loading — this is the dominant band-broadening
mechanism and it sets your resolution limit.

### The paradox worth internalizing

$D_{eff} \propto 1/D$. **Faster-diffusing molecules disperse LESS.** Increase $D$
and molecules cross streamlines quickly, sampling all velocities and averaging
out the shear. Decrease $D$ and they stay stuck on their streamline, experiencing
the full velocity spread.

So a small ion spreads *less* along the channel than a protein does, even though
it diffuses faster across it. This runs against intuition and is a genuinely
useful thing to know when designing a separation.

### Validity

The formula assumes enough time for cross-stream equilibration:
$t \gg h^2/D$. Before that you are in the pre-asymptotic regime and dispersion
is *less* than the formula predicts. Check:

$$t_{transit} = \frac{L}{\bar u} \gg \frac{h^2}{D}$$

`taylor_dispersion.m` verifies both the formula and this validity condition
numerically, by tracking the variance of a plug and comparing to $2D_{eff}t$.

---

## 9.6 Making it mix: the three strategies

### 1. Just make the channel longer

Works. Costs area. $L \approx 0.3\,Pe\,w$.

### 2. Lamination — cut the diffusion distance

Split and recombine to make $n$ thinner stripes. Diffusion time scales as
(stripe width)², so $n$ stripes mix $n^2$ times faster:

$$t_{mix} \sim \frac{(w/n)^2}{D}$$

A 4-stripe laminator mixes **16×** faster. This is why interdigitated inlet
manifolds exist, and it is the highest-value idea in passive mixing.

### 3. Chaotic advection — stretch and fold

The clever one. Stokes flow is reversible, but *chaotic* Stokes flow still
stretches material lines **exponentially**:

$$\ell(t) \sim \ell_0 e^{\lambda t}$$

where $\lambda$ is the Lyapunov exponent. Stripe width shrinks exponentially, so
mixing time grows only **logarithmically** with $Pe$:

$$t_{mix} \sim \frac{1}{\lambda}\ln(Pe)$$

versus $t_{mix}\sim Pe$ for pure diffusion. At $Pe = 10^4$ that is the
difference between a metre of channel and a millimetre.

The canonical implementation is the **staggered herringbone mixer** (Stroock et
al., *Science* 2002): asymmetric grooves in the channel floor drive a transverse
helical flow, and alternating the groove asymmetry every half-cycle makes the
flow chaotic. Crucially, it is a **3D** effect — the transverse flow lives in the
channel cross-section — so **a 2D simulation cannot reproduce it.** A 2D
serpentine gets you a mild stretching effect and nothing more; do not conclude
from a 2D run that herringbones don't help.

**What 2D *can* show you** is the serpentine's genuine 2D effect: at bends, the
inside of the turn moves faster, which shears the interface and increases
interfacial length. `mixing_study.m` measures exactly this, comparing the
serpentine against a straight channel of the same length.

---

## 9.7 The workflow

```matlab
% 1. Flow field (once)
serpentine                      % saves serpentine_flow.mat

% 2. Transport the scalar on it
mixing_study                    % loads it, runs both geometries
```

`mixing_study.m` reports, for each geometry:

- $M(x)$ down the channel
- interfacial length $\int|\nabla c|$
- distance to $M < 0.05$
- **measured $D_{num}/D_{phys}$** for the grid used
- a grid-refinement check

That fourth item is the one that makes the result trustworthy. The others are
just numbers until you know the mesh isn't producing them.

---

## Exercises

**9.1** Run `mixing_study`. Compare the mixing length of the serpentine to the
straight channel. Is the improvement large? *(It should be modest — this is the
honest 2D result, and it is why real mixers use 3D geometries.)*

**9.2** Change the scheme from `'vanleer'` to `'upwind'` and rerun. How much
shorter does the mixing length appear? That difference is entirely numerical
diffusion, and it is the size of the error you would have published.

**9.3** Run `taylor_dispersion`. Confirm the measured $D_{eff}$ matches
$D(1+2Pe_h^2/105)$. Then reduce the channel length until $t_{transit} < h^2/D$
and observe the pre-asymptotic regime, where the formula overpredicts.

**9.4 (the paradox)** Compute $D_{eff}$ for a small ion ($D = 2\times10^{-9}$)
and a protein ($D = 4\times10^{-11}$) in the same channel. Which plug spreads
more in absolute terms? Verify numerically.

**9.5** Implement lamination: modify the inlet condition in `mixing_study` to
create 2, 4, and 8 alternating stripes. Confirm mixing length scales as $1/n^2$.

**9.6** Compute the Lyapunov exponent of your serpentine flow by advecting two
initially-adjacent passive tracers and tracking their separation. Is it positive?
*(For a 2D steady flow it should be essentially zero — steady 2D flow cannot be
chaotic, because streamlines are level sets of the streamfunction and particles
cannot cross them. This is a theorem, and confirming it numerically explains
precisely why chaotic micromixers must be either 3D or time-dependent.)*

**9.7 (design)** You must mix a protein solution in a 50 µm channel at 2 mm/s.
Budget: 10 mm of channel. Using the numbers above, decide whether a plain
channel, a serpentine, a 4-way laminator, or a herringbone is required. Justify
with calculations, and state which parts your 2D code can verify and which it
cannot.

---

**Next:** [Tutorial 10 — Surface tension and capillarity](10-surface-tension-capillarity.md).
Part III turns to the second phenomenon you asked for: interfaces.
