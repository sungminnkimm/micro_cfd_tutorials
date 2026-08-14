function opt = ns2_defaults()
%NS2_DEFAULTS  Default options for NS2PHASE (two-phase solver).
%
%   Defaults are a WATER-IN-OIL system, because a density ratio near 1 is much
%   kinder numerically than water/air's 1000:1 and is what droplet
%   microfluidics actually uses.  Start here; move to high density ratios only
%   once everything else is validated.
%
%   Phase 1 (phi > 0) is the DISPERSED phase (the droplet).
%   Phase 2 (phi < 0) is the CONTINUOUS phase.
%
%   FIELDS beyond NS_DEFAULTS:
%     rho1, mu1     dispersed-phase properties      998, 1.0e-3
%     rho2, mu2     continuous-phase properties     900, 1.0e-2
%     sigma         interfacial tension [N/m]       5e-3
%     phi0          initial level set (REQUIRED)    []
%     epsFac        interface half-width / dx       1.5
%     reinitEvery   reinitialize every N steps      10
%     reinitIters   pseudo-time iterations          5
%     Uref          reference speed for the CFL dt  1e-3
%
%   opt.solid, opt.uMask, opt.vMask: set with SET_SOLID before calling.
%
%   See also NS2PHASE, LEVEL_SET, SET_SOLID.

opt.nx = 64;    opt.ny = 64;
opt.Lx = 200e-6;  opt.Ly = 200e-6;

opt.rho1 = 998;   opt.mu1 = 1.0e-3;    % water droplet
opt.rho2 = 900;   opt.mu2 = 1.0e-2;    % oil continuous phase
opt.sigma = 5e-3;                       % water/oil + surfactant

opt.gx = 0;  opt.gy = 0;
opt.Uref = 1e-3;

opt.tEnd  = 1e-3;
opt.cfl   = 0.4;
opt.dtMax = Inf;
opt.report = 100;

opt.phi0 = [];
opt.u0 = [];  opt.v0 = [];
opt.epsFac = 1.5;
opt.reinitEvery = 20;
opt.reinitIters = 3;

% conserveMass: after each reinitialization, shift phi to restore the INITIAL
% enclosed area.  Essential for a closed domain (a static droplet otherwise
% dissolves).  Set FALSE when the domain has inflow/outflow and the dispersed
% volume is genuinely supposed to change -- e.g. droplet generation, where the
% whole point is that new fluid keeps arriving.
opt.conserveMass = true;

opt.bc.left   = struct('type','wall');
opt.bc.right  = struct('type','wall');
opt.bc.bottom = struct('type','wall');
opt.bc.top    = struct('type','wall');

opt.solid = [];  opt.uMask = [];  opt.vMask = [];

% Refuse to start a run needing more than this many steps.  A capillary- or
% viscous-limited two-phase case can silently demand millions: the first draft
% of T_JUNCTION asked for 5.3 million (~90 min) and gave no warning until it
% failed to finish.  Raise this deliberately if you really want a long run.
opt.maxSteps = 3e5;
end
