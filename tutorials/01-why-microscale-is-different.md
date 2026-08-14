# Tutorial 1 — Why Microscale Is Different

**Goal:** Understand which terms in the Navier–Stokes equations you are allowed
to drop at the microscale, and why. Build a scaling calculator you will use for
the rest of the course.

**Code:** `code/ch01/scaling.m`, `code/ch01/scaling_demo.m`, `code/ch01/scale_sweep.m`

---

## 1.1 The one idea

Shrinking a device does not scale all physics equally. Every physical effect has
its own power of $L$:

| Effect | Scales as | At $L = 100\ \mu$m vs 1 m |
|--------|-----------|----------------------------|
| Inertia | $L^3$ (mass × accel per unit … ) | $10^{-12}$ |
| Viscous force | $L^1$ | $10^{-4}$ |
| Gravity / buoyancy | $L^3$ | $10^{-12}$ |
| Surface tension | $L^1$ | $10^{-4}$ |
| Diffusion time | $L^2$ | $10^{-8}$ |
| Surface area / volume | $L^{-1}$ | $10^{4}$ |

Inertia and gravity fall off as $L^3$. Viscosity and surface tension fall off as
$L^1$. So as you shrink, *viscosity and surface tension win*, and they win by
four orders of magnitude by the time you reach 100 µm.

Everything else in this tutorial is bookkeeping on that one sentence.

---

## 1.2 Non-dimensionalizing Navier–Stokes

Start from incompressible Navier–Stokes:

$$
\rho\left(\frac{\partial \mathbf{u}}{\partial t} + \mathbf{u}\cdot\nabla\mathbf{u}\right)
= -\nabla p + \mu \nabla^2 \mathbf{u} + \rho \mathbf{g}
$$

Introduce scales: length $L$, velocity $U$, time $T = L/U$, pressure $P$.
Write $\mathbf{u} = U\mathbf{u}^*$, $\nabla = \frac{1}{L}\nabla^*$, $t = \frac{L}{U}t^*$.
Substituting and dividing through by $\mu U / L^2$ (the *viscous* scale — this is
the right choice at the microscale, and choosing it is half the insight):

$$
\underbrace{\frac{\rho U L}{\mu}}_{Re}\left(\frac{\partial \mathbf{u}^*}{\partial t^*} + \mathbf{u}^*\cdot\nabla^*\mathbf{u}^*\right)
= -\nabla^* p^* + \nabla^{*2}\mathbf{u}^* + \underbrace{\frac{\rho g L^2}{\mu U}}_{\text{Stokes/gravity}}\hat{\mathbf{g}}
$$

with $P = \mu U / L$ — note that microscale pressure scales *viscously*, not as
$\rho U^2$. The dynamic-pressure scale $\rho U^2$ is meaningless here; if you
size a pump using it you will be wrong by a factor of $1/Re$, i.e. by ten.
Viscous scaling is also why pressure drop is so brutally sensitive to channel
size: $\Delta p/\Delta x = 12\mu U/h^2 \propto h^{-2}$ at fixed velocity, and
$\propto h^{-4}$ at fixed *flow rate*. Halving a channel's width at fixed flow
rate costs you 16× the pressure.

**The Reynolds number is now a coefficient multiplying the entire left-hand side.**
If $Re \ll 1$, the left-hand side is negligible and you get the **Stokes equations**:

$$
0 = -\nabla p + \mu\nabla^2\mathbf{u}, \qquad \nabla\cdot\mathbf{u} = 0
$$

This is a *linear*, *time-independent*, *instantaneous* system. Three enormous
consequences:

1. **Linearity.** Superposition works. Double the pressure drop, double the flow.
2. **No history.** The flow field depends only on the *current* boundary
   conditions. Stop pushing and the fluid stops immediately (in ~$\rho L^2/\mu$
   ≈ 10 µs for water at 100 µm).
3. **Reversibility.** Run the boundary motion backwards and the fluid retraces
   its path exactly. This is the famous dye-in-glycerin demonstration, and it is
   *why microchannels refuse to mix*.

Consequence 3 is the central pain of microfluidics and the subject of Tutorial 9.

---

## 1.3 The dimensionless numbers you actually need

### Reynolds number — inertia vs viscosity

$$Re = \frac{\rho U L}{\mu}$$

- $Re \lesssim 1$: Stokes flow. Steady, reversible, linear.
- $1 \lesssim Re \lesssim 100$: inertia matters but flow is laminar and steady.
  Inertial microfluidics lives here (particle focusing, Dean flow).
- $Re \gtrsim 2000$: turbulence. You will not see this in a microchannel without
  trying very hard.

**Numerically:** low $Re$ means the convective term is weak, which means
central differencing is stable and you don't need upwinding for momentum. It
also means the diffusive timestep limit dominates — see Tutorial 4.

### Péclet number — advection vs diffusion (of a scalar)

$$Pe = \frac{U L}{D}$$

Here $D$ is *mass* diffusivity, not momentum diffusivity. This is the number
that surprises people: $Re$ is tiny but $Pe$ is huge, because for liquids
$\nu/D = Sc \sim 10^3$. Water at 1 mm/s in a 100 µm channel:

$$Re = 0.1, \qquad Pe = \frac{10^{-3}\times 10^{-4}}{5\times10^{-10}} = 200$$

So momentum diffuses across the channel instantly, but *molecules do not*. Two
streams flowing side by side stay side by side. The diffusive mixing length is

$$L_{mix} \sim Pe \cdot w = 200 \times 100\,\mu\text{m} = 2\ \text{cm}$$

of channel — for a 100 µm channel. That's the whole reason mixer geometries exist.

**Numerically:** high $Pe$ means the scalar transport equation is advection-dominated,
which means you *do* need upwinding or a high-resolution scheme, and you will
fight numerical diffusion. Tutorial 5.

### Capillary number — viscous vs surface tension

$$Ca = \frac{\mu U}{\sigma}$$

- $Ca \ll 10^{-3}$: interfaces are governed by surface tension. Droplets are
  spheres/plugs, menisci are circular arcs. Squeezing regime in a T-junction.
- $Ca \sim 10^{-2}$: viscous shear starts deforming interfaces. Dripping.
- $Ca \gtrsim 10^{-1}$: shear wins, interfaces stretch into jets/threads.

Water at 1 mm/s: $Ca = 1.4\times10^{-5}$. Surface tension dominates absolutely.

**Numerically:** low $Ca$ is *hard*. Strong surface tension means a stiff
capillary wave timestep restriction and severe spurious currents. Tutorial 11.

### Bond number — gravity vs surface tension

$$Bo = \frac{\rho g L^2}{\sigma}$$

At 100 µm: $Bo = \frac{998 \times 9.81 \times 10^{-8}}{0.072} = 1.4\times10^{-3}$.

**Gravity is irrelevant.** Droplets don't sag, don't settle appreciably, don't
care which way is up. You can almost always drop the gravity term. The crossover
is the **capillary length** $\ell_c = \sqrt{\sigma/\rho g} \approx 2.7$ mm for
water — anything smaller than that is surface-tension-dominated.

### Weber number — inertia vs surface tension

$$We = \frac{\rho U^2 L}{\sigma} = Re \cdot Ca$$

Tiny at microscale. Droplets don't splash.

### Knudsen number — continuum validity

$$Kn = \frac{\lambda}{L}$$

$\lambda$ = molecular mean free path (~68 nm for air at STP, ~0.3 nm for liquid
water). Continuum Navier–Stokes with no-slip is valid for $Kn < 10^{-3}$.

- **Liquids:** valid down to ~10 nm. You are always fine.
- **Gases:** at $L = 1\ \mu$m, $Kn = 0.068$ — you are in the slip regime and need
  slip boundary conditions.

For this course (liquid microfluidics) $Kn$ never matters. It is here so you
know when to stop trusting the code.

---

## 1.4 Build the calculator

`code/ch01/scaling.m` — a function that takes a flow condition and returns every
dimensionless group, plus the derived time and length scales you'll need to size
your grid and timestep.

```matlab
function s = scaling(varargin)
%SCALING  Dimensionless groups and derived scales for a microscale flow.
%
%   s = scaling('L',100e-6, 'U',1e-3)                 % water defaults
%   s = scaling('L',50e-6,  'U',5e-3, 'D',1e-10)
%   s = scaling(...,'rho',1000,'mu',1e-3,'sigma',0.072,'g',9.81)
%
%   Returns a struct with fields Re, Pe, Ca, Bo, We, Sc, plus:
%     t_visc   momentum diffusion time across L       (rho*L^2/mu)
%     t_diff   mass diffusion time across L           (L^2/D)
%     t_adv    advection time across L                (L/U)
%     t_cap    capillary wave time on scale L         (sqrt(rho*L^3/sigma))
%     L_mix    channel length for diffusive mixing    (Pe*L)
%     l_cap    capillary length                       (sqrt(sigma/(rho*g)))
%     dp_pois  pressure gradient for this U in a slot (12*mu*U/L^2)

p = inputParser;
p.addParameter('L',     100e-6);   % m      characteristic length
p.addParameter('U',     1e-3);     % m/s    characteristic velocity
p.addParameter('rho',   998);      % kg/m^3
p.addParameter('mu',    1.0e-3);   % Pa.s
p.addParameter('sigma', 0.072);    % N/m
p.addParameter('D',     5e-10);    % m^2/s  small-molecule in water
p.addParameter('g',     9.81);     % m/s^2
p.parse(varargin{:});
q = p.Results;

L = q.L;  U = q.U;  rho = q.rho;  mu = q.mu;  sig = q.sigma;  D = q.D;  g = q.g;
nu = mu/rho;

s        = q;
s.nu     = nu;

% --- dimensionless groups ------------------------------------------------
s.Re = rho*U*L/mu;          % inertia   / viscous
s.Pe = U*L/D;               % advection / mass diffusion
s.Ca = mu*U/sig;            % viscous   / surface tension
s.Bo = rho*g*L^2/sig;       % gravity   / surface tension
s.We = rho*U^2*L/sig;       % inertia   / surface tension  (= Re*Ca)
s.Sc = nu/D;                % momentum diffusivity / mass diffusivity

% --- time scales ---------------------------------------------------------
s.t_visc = L^2/nu;              % momentum diffuses across L
s.t_diff = L^2/D;               % molecules diffuse across L
s.t_adv  = L/U;                 % fluid crosses L
s.t_cap  = sqrt(rho*L^3/sig);   % capillary wave period on scale L

% --- length scales -------------------------------------------------------
s.L_mix = s.Pe*L;                    % channel length to mix by diffusion
s.l_cap = sqrt(sig/(rho*g));         % capillary length
s.d_ent = 0.06*s.Re*L;               % hydrodynamic entrance length (laminar)

% --- engineering ---------------------------------------------------------
s.dp_pois = 12*mu*U/L^2;             % dP/dx for mean U in a slot of height L
end
```

And a reporter, `code/ch01/scaling_demo.m`:

```matlab
%SCALING_DEMO  Print the dimensionless picture for a typical microchannel.
if ~exist('QUIET','var'), QUIET = false; end

s = scaling('L',100e-6,'U',1e-3);

fprintf('\n--- Microchannel: L = %g um, U = %g mm/s (water) ---\n', ...
        s.L*1e6, s.U*1e3);
fprintf('  Re = %10.3g   %s\n', s.Re, ch01_verdict(s.Re, [1 100], ...
        'Stokes flow: drop inertia', 'inertia matters', 'transitional'));
fprintf('  Pe = %10.3g   %s\n', s.Pe, ch01_verdict(s.Pe, [1 100], ...
        'diffusion-dominated', 'advection-dominated: mixing is HARD', 'comparable'));
fprintf('  Ca = %10.3g   %s\n', s.Ca, ch01_verdict(s.Ca, [1e-3 1e-1], ...
        'surface tension dominates interfaces', 'shear deforms interfaces', 'transitional'));
fprintf('  Bo = %10.3g   %s\n', s.Bo, ch01_verdict(s.Bo, [1e-2 1], ...
        'gravity negligible', 'gravity matters', 'comparable'));
fprintf('  We = %10.3g\n', s.We);
fprintf('  Sc = %10.3g   (momentum diffuses %.0fx faster than mass)\n', s.Sc, s.Sc);

fprintf('\n  Time scales:\n');
fprintf('    momentum diffusion across L : %10.3g s\n', s.t_visc);
fprintf('    advection across L          : %10.3g s\n', s.t_adv);
fprintf('    capillary wave              : %10.3g s\n', s.t_cap);
fprintf('    MASS diffusion across L     : %10.3g s   <-- the slow one\n', s.t_diff);

fprintf('\n  Design numbers:\n');
fprintf('    channel length to mix by diffusion : %8.3g mm\n', s.L_mix*1e3);
fprintf('    capillary length (water)           : %8.3g mm\n', s.l_cap*1e3);
fprintf('    entrance length                    : %8.3g um\n', s.d_ent*1e6);
fprintf('    dP/dx to drive this flow           : %8.3g bar/cm\n', s.dp_pois*1e-5*1e-2);
fprintf('\n');
```

with the small helper in its own file, `code/ch01/ch01_verdict.m`:

```matlab
function txt = ch01_verdict(val, lims, lowTxt, highTxt, midTxt)
if     val < lims(1), txt = ['<-- ' lowTxt];
elseif val > lims(2), txt = ['<-- ' highTxt];
else,                 txt = ['<-- ' midTxt];
end
end
```

Run it:

```matlab
scaling_demo
```

You should see something close to:

```
--- Microchannel: L = 100 um, U = 1 mm/s (water) ---
  Re =     0.0998   <-- Stokes flow: drop inertia
  Pe =        200   <-- advection-dominated: mixing is HARD
  Ca =   1.39e-05   <-- surface tension dominates interfaces
  Bo =    0.00136   <-- gravity negligible
  We =   1.39e-06
  Sc =      2e+03   (momentum diffuses 2004x faster than mass)

  Time scales:
    momentum diffusion across L :    0.00998 s
    advection across L          :        0.1 s
    capillary wave              :   0.000118 s
    MASS diffusion across L     :         20 s   <-- the slow one

  Design numbers:
    channel length to mix by diffusion :       20 mm
    capillary length (water)           :     2.71 mm
    entrance length                    :    0.599 um
    dP/dx to drive this flow           :  0.00012 bar/cm
```

Two of those design numbers deserve a second look. The **entrance length is 0.6
µm** — six *tenths* of a micron. The flow becomes fully developed essentially at
the inlet, so you can use the analytic Poiseuille profile as an inlet boundary
condition without apology (Tutorial 8 does exactly this). And the **mixing length
is 20 mm** in a device that is probably 20 mm long in total. Diffusive mixing
will consume your entire chip.

**Sit with that for a moment.** $Re = 0.1$ and $Pe = 200$ describe *the same
flow*. Momentum is diffusion-dominated; mass is advection-dominated. That
mismatch — a factor of $Sc = 2000$ — is the defining feature of liquid
microfluidics, and it is why your solver will need two different numerical
treatments for the two equations.

---

## 1.5 The sweep: seeing the crossovers

`code/ch01/scale_sweep.m` plots each group against channel size, so you can *see*
where each physical effect switches off.

```matlab
%SCALE_SWEEP  How the dimensionless groups vary with device size.
if ~exist('QUIET','var'), QUIET = false; end

L = logspace(-7, -2, 200);      % 100 nm .. 1 cm
U = 1e-3;                       % fixed 1 mm/s

Re = zeros(size(L)); Pe = Re; Ca = Re; Bo = Re;
for k = 1:numel(L)
    s = scaling('L',L(k),'U',U);
    Re(k)=s.Re; Pe(k)=s.Pe; Ca(k)=s.Ca; Bo(k)=s.Bo;
end

% crossover sizes (where each group = 1)
L_Re = interp1(log(Re), L, 0, 'linear', NaN);
L_Pe = interp1(log(Pe), L, 0, 'linear', NaN);
L_Bo = interp1(log(Bo), L, 0, 'linear', NaN);

fprintf('At U = %g mm/s:\n', U*1e3);
fprintf('  Re = 1 at L = %8.3g mm  (above: inertia matters)\n', L_Re*1e3);
fprintf('  Pe = 1 at L = %8.3g um  (above: advection beats diffusion)\n', L_Pe*1e6);
fprintf('  Bo = 1 at L = %8.3g mm  (above: gravity beats surface tension)\n', L_Bo*1e3);

if ~QUIET
    figure('Color','w','Position',[100 100 760 460]);
    loglog(L*1e6, Re,'LineWidth',2); hold on
    loglog(L*1e6, Pe,'LineWidth',2);
    loglog(L*1e6, Ca,'LineWidth',2);
    loglog(L*1e6, Bo,'LineWidth',2);
    yline(1,'k--','= 1','LineWidth',1.2);
    xregion = [10 500];
    patch([xregion fliplr(xregion)], [1e-10 1e-10 1e8 1e8], ...
          [0.2 0.5 0.9],'FaceAlpha',0.07,'EdgeColor','none');
    text(70, 1e6, 'typical microfluidics','Color',[0.2 0.4 0.7]);
    grid on; xlabel('channel size L  [\mum]'); ylabel('dimensionless group');
    legend({'Re','Pe','Ca','Bo'},'Location','northwest');
    title(sprintf('Scaling at U = %g mm/s (water)', U*1e3));
    ylim([1e-10 1e8]);
end
```

Read the plot. In the shaded microfluidic band: $Pe$ is the only group above 1.
That is the whole story.

---

## 1.6 What this buys you numerically

This is the part that matters for the code you're about to write.

| Physical fact | Numerical consequence |
|---|---|
| $Re \ll 1$ | Convective term is weak → central differencing on momentum is fine, no upwinding needed. But the CFL limit is *not* what constrains you. |
| $Re \ll 1$ | The **viscous** timestep limit $\Delta t < \frac{1}{2\nu}\left(\frac{1}{\Delta x^2}+\frac{1}{\Delta y^2}\right)^{-1}$ dominates. Explicit diffusion gets expensive fast. Tutorial 4. |
| $Pe \gg 1$ | Scalar transport *is* advection-dominated → you need upwind/TVD schemes, and numerical diffusion will masquerade as physical mixing. Tutorial 5 and 9. |
| $Bo \ll 1$ | Drop gravity. One less term, one less source of bugs. |
| $Ca \ll 1$ | Surface tension is stiff → capillary timestep $\Delta t < \sqrt{\frac{\bar\rho \Delta x^3}{2\pi\sigma}}$. This will be your smallest timestep in Tutorial 11. |
| Steady Stokes | You can solve for the flow field *once* and freeze it, then transport scalars on it for as long as you like. Huge savings in Tutorial 9. |

That last row is worth stating plainly: **in a rigid microchannel with steady
boundary conditions, the velocity field is constant in time.** You do not need to
re-solve Navier–Stokes every step to study mixing. Solve once, save `u` and `v`,
then run the scalar equation. Tutorial 9 does exactly this.

---

## Exercises

**1.1** Use `scaling` to find the velocity at which $Re = 1$ in a 50 µm channel.
Is that a velocity you could achieve with a syringe pump? (Typical range: 0.1
µL/min to 10 mL/min.)

**1.2** A protein has $D \approx 4\times10^{-11}$ m²/s. Recompute $Pe$ and $L_{mix}$
for the standard microchannel. How long must the channel be? Now do the same for
a 1 µm bacterium ($D \approx 4\times10^{-13}$ m²/s). What does this tell you
about separating cells by diffusion?

**1.3** You want to make 50 µm water droplets in oil ($\mu_{oil} = 10$ mPa·s,
$\sigma = 5$ mN/m) at 1 mm/s. Compute $Ca$ using the *oil* viscosity. Which
T-junction regime does Tutorial 12's map put you in?

**1.4** Modify `scale_sweep.m` to sweep *velocity* at fixed $L = 100$ µm instead.
At what velocity does $Ca$ reach $10^{-2}$? Is that flow still laminar?

**1.5 (the one that matters)** Take a device you actually care about — real
dimensions, real flow rate, real fluids. Run `scaling` on it. Write down, in one
sentence each: which terms of Navier–Stokes you can drop, whether you need a
two-phase model, and how long the channel must be to mix. Keep that note. You
will check your Tutorial 8 and 9 results against it.

---

**Next:** [Tutorial 2 — Finite differences and truncation error](02-finite-differences.md).
Before you can solve Navier–Stokes you need to be able to differentiate a
function on a grid and *prove* you did it correctly.
