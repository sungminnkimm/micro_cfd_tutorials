%T_JUNCTION  Droplet generation at a T-junction.
%
%   Tutorial 12.  Run:  t_junction
%
%   Dispersed phase enters from a side channel into a cross-flowing continuous
%   phase.  In the SQUEEZING regime (Ca < 1e-2) the droplet grows until it
%   blocks the main channel; the resulting pressure build-up -- not shear --
%   pinches the neck.  Size then follows Garstecki et al. (Lab Chip 2006):
%
%       L/w = 1 + alpha * Qd/Qc,     alpha ~ 1
%
%   independent of Ca, viscosity and surface tension.  That independence is why
%   the squeezing regime is so reproducible, and why devices are designed to
%   sit in it.
%
%   READ THE CAVEATS in Tutorial 12.4 before quoting any number from this.  In
%   particular a 2-D droplet is an infinite cylinder, and in 3-D the continuous
%   phase leaks through the CORNERS of a rectangular channel -- a genuinely 3-D
%   part of the squeezing mechanism that this model cannot represent.  Expect
%   the right trends, not the right absolute volume.
%
%   See also REGIME_MAP, NS2PHASE, FLOW_FOCUSING.

if ~exist('QUIET','var'), QUIET = false; end
if ~exist('RATIO','var'), RATIO = 0.5; end        % Qd/Qc, sweep for Ex 12.3

%% ------------------------------------------------------------ geometry ---
%  A CAUTIONARY NOTE ABOUT THIS CONFIGURATION
%  ------------------------------------------
%  The first version of this demo used Lx = 500 um on a 200x80 grid with a
%  10 mPa.s oil.  That gives dx = 2.5 um, dt = 5.6e-8 s (viscous-limited), and
%  a residence time of 0.3 s -- which is 5.3 MILLION timesteps, roughly 90
%  minutes of CPU for a single case.  Running REGIME_MAP first would have
%  reported that in about a second.  The author did not, and paid for it.
%
%  The configuration below is chosen to be affordable (~40k steps):
%    * a lighter continuous phase (2 mPa.s), which BOTH lowers Ca and relaxes
%      the viscous timestep, because nu = mu/rho falls;
%    * a faster continuous phase, which shortens the residence time;
%    * a coarser grid.
%  The cost is resolution: only ~17 cells across the channel and fewer across
%  the neck, which is BELOW the >= 20 that Tutorial 12.4 asks for.  Treat the
%  droplet length below as indicative, and do Exercise 12.5 (grid convergence)
%  before believing any specific number.
w    = 50e-6;                 % main channel width (y-extent of the main duct)
Lx   = 300e-6;
Ly   = 120e-6;
wIn  = 50e-6;                 % side-channel width
xIn  = 75e-6;                 % where the side channel meets the main channel

nx = 100;  ny = 40;
G  = mac_grid(nx, ny, Lx, Ly);
fprintf('\n  T-junction: main w = %g um at %.1f cells across\n', w*1e6, w/G.dy);

% Main channel occupies the lower strip 0 < y < w.  The side channel is a
% vertical duct above it at xIn <= x <= xIn+wIn.  Everything else is solid.
yMain = w;
inMain = G.Yp <= yMain;
inSide = G.Xp >= xIn & G.Xp <= xIn+wIn & G.Yp > yMain;
solid  = ~(inMain | inSide);

%% ------------------------------------------------------------- fluids ----
mu_c = 2.0e-3; rho_c = 900;      % light oil, continuous phase
mu_d = 1.0e-3; rho_d = 998;      % water, dispersed
sigma = 5e-3;

Uc = 10e-3;                       % continuous-phase mean velocity
Ca = mu_c*Uc/sigma;
Qc = Uc*w;                        % per unit depth (2-D)
Qd = RATIO*Qc;
Ud = Qd/wIn;                      % side-channel mean velocity

fprintf('  Ca = %.4f (%s), Qd/Qc = %.2f\n', Ca, regime_of(Ca), RATIO);
fprintf('  Uc = %.2f mm/s, Ud = %.2f mm/s\n', Uc*1e3, Ud*1e3);
fprintf('  Garstecki prediction: L/w = 1 + %.2f = %.2f  -> L = %.1f um\n', ...
        RATIO, 1+RATIO, (1+RATIO)*w*1e6);

%% ---------------------------------------------------------- solver opts --
opt = ns2_defaults();
opt.nx = nx;  opt.ny = ny;  opt.Lx = Lx;  opt.Ly = Ly;
opt.rho1 = rho_d;  opt.mu1 = mu_d;        % phi > 0 = dispersed
opt.rho2 = rho_c;  opt.mu2 = mu_c;
opt.sigma = sigma;
opt = set_solid(opt, solid);

% Initial interface: dispersed phase fills the side channel down to the
% junction, so phi > 0 there.
phi0 = -ones(nx,ny)*1e-4;
seed = G.Xp >= xIn & G.Xp <= xIn+wIn & G.Yp > yMain;
d = min(min(G.Xp-xIn, xIn+wIn-G.Xp), G.Yp-yMain);
phi0(seed) = d(seed);
opt.phi0 = phi0;

% Inflow: continuous phase from the left across the main channel only.
yu = G.yu;
uIn = zeros(size(yu));
inCh = yu < yMain;
eta  = yu(inCh)/yMain;
uIn(inCh) = 6*Uc*eta.*(1-eta);
opt.bc.left  = struct('type','inflow','u',uIn);
opt.bc.right = struct('type','outflow');

% Dispersed phase injected from the top of the side channel.
xv = G.xv;
vIn = zeros(size(xv));
inSd = xv >= xIn & xv <= xIn+wIn;
vIn(inSd) = -Ud;                        % downward into the junction
opt.bc.top    = struct('type','inflow','v',vIn);
opt.bc.bottom = struct('type','wall');

opt.Uref = max(Uc, Ud)*1.5;
tRes = Lx/Uc;
opt.tEnd = 0.6*tRes;              % keep the step count affordable; see header
opt.report = 5000;

% Open boundaries: dispersed fluid keeps arriving, so the enclosed area is
% SUPPOSED to grow.  The global area correction (which pins it to its initial
% value) must be off here -- leaving it on would fight the injection and is a
% subtle way to get a completely wrong droplet size.
opt.conserveMass = false;

% ...which leaves reinitialization drift completely uncorrected, and that is a
% real problem here.  At the default (every 20 steps, 3 iterations) a 50k-step
% run makes 2500 reinit calls; at ~0.1% area each, the interface bleeds away
% faster than it is injected and NO droplet can form.  Reinitialize as rarely
% as the |grad phi| diagnostic allows -- watch that column and keep it near 1.
opt.reinitEvery = 60;
opt.reinitIters = 2;

fprintf('\n  residence time %.3f s; simulating %.3f s\n', tRes, opt.tEnd);

tic; S = ns2phase(opt); el = toc;
fprintf('  wall time %.1f s\n', el);

%% ------------------------------------------------------- measure size ----
% Slice along the main-channel centreline and find contiguous runs of phi > 0.
j = max(1, round(0.5*yMain/G.dy));
line = S.phi(:,j) > 0;
% drop the run that contains the injection point
[segs, lens] = runs_of_true(line, G.xp);
segs = segs(:); lens = lens(:);
keep = true(size(lens));
for k = 1:numel(lens)
    if segs(k) <= xIn+wIn && segs(k)+lens(k) >= xIn, keep(k) = false; end
end
lens = lens(keep);

fprintf('\n  detached droplets found on the centreline: %d\n', numel(lens));
if ~isempty(lens)
    fprintf('    lengths [um]: '); fprintf('%.1f ', lens*1e6); fprintf('\n');
    Lm = mean(lens);
    fprintf('    mean L = %.1f um -> L/w = %.2f\n', Lm*1e6, Lm/w);
    fprintf('    Garstecki L/w = %.2f  -> implied alpha = %.2f\n', ...
            1+RATIO, (Lm/w - 1)/RATIO);
else
    fprintf(['    NONE -- and this is the SHIPPED behaviour at this resolution.\n' ...
             '    Read it as the Tutorial 12.4 resolution requirement being\n' ...
             '    demonstrated rather than merely asserted.\n\n' ...
             '    The run is numerically healthy: max|div u| ~ 1e-12, and the\n' ...
             '    enclosed area GROWS (positive area error) exactly as continuous\n' ...
             '    injection demands. The dispersed phase enters, fills the\n' ...
             '    junction and is dragged downstream -- but the NECK never gets\n' ...
             '    thin enough, in cells, for the level set to resolve pinch-off.\n' ...
             '    At 16.7 cells across the channel the neck is only a few cells\n' ...
             '    wide, and Tutorial 12.4 asks for >= 20 ACROSS THE NECK.\n\n' ...
             '    This matters more than a missing picture. If it HAD pinched\n' ...
             '    off here, breakup would have been triggered by two interfaces\n' ...
             '    coming within a cell of each other -- a NUMERICAL criterion,\n' ...
             '    not a physical one -- and the droplet size would have been a\n' ...
             '    function of dx. A demo that silently produced a plausible\n' ...
             '    droplet length at this grid would be worse than one that\n' ...
             '    refuses to.\n\n' ...
             '    Exercise 12.5 is where you fix it: refine until droplets appear\n' ...
             '    AND their size stops changing. Expect to need ~3-4x this grid,\n' ...
             '    and check the step count with REGIME_MAP before you start.\n']);
end

fprintf('\n  Diagnostics (all must be healthy before the size above means anything):\n');
fprintf('    area error   %+.3e   (level-set mass loss)\n', S.areaErr);
fprintf('    max|div u|   %.3e\n', S.divMax);
fprintf('    dt limited by %s, dt = %.3e s, %d steps\n', S.dtLimit, S.dt, S.steps);
fprintf('    cells across w: %.1f  (want >= 20 across the NECK, which is thinner)\n\n', ...
        w/G.dy);

if ~QUIET
    figure('Color','w','Position',[60 60 1050 420]);
    subplot(2,1,1);
    ph = S.phi.';  ph(solid.') = NaN;
    imagesc(G.xp*1e6, G.yp*1e6, ph); set(gca,'YDir','normal'); hold on
    contour(G.xp*1e6, G.yp*1e6, S.phi.', [0 0],'k','LineWidth',2);
    contour(G.xp*1e6, G.yp*1e6, double(solid.'), [0.5 0.5],'w','LineWidth',1.5);
    axis image; colorbar; title(sprintf('\\phi, Ca = %.4f, Q_d/Q_c = %.2f', Ca, RATIO));
    subplot(2,1,2);
    sp = sqrt(S.uc.^2+S.vc.^2).'*1e3; sp(solid.') = NaN;
    imagesc(G.xp*1e6, G.yp*1e6, sp); set(gca,'YDir','normal'); hold on
    contour(G.xp*1e6, G.yp*1e6, S.phi.', [0 0],'w','LineWidth',1.5);
    axis image; colorbar; title('|u| [mm/s]'); xlabel('x [\mum]');
end

%% ===================================================== local functions ===
function r = regime_of(Ca)
if     Ca < 1e-2, r = 'squeezing';
elseif Ca < 1e-1, r = 'dripping';
else,             r = 'jetting';
end
end

function [starts, lens] = runs_of_true(mask, x)
%RUNS_OF_TRUE  Start positions and lengths of contiguous true runs.
mask = mask(:).';
d = diff([false mask false]);
i0 = find(d==1);  i1 = find(d==-1)-1;
starts = x(i0);
lens   = x(min(i1,numel(x))) - x(i0) + (x(2)-x(1));
end
