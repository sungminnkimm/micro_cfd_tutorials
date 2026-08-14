function opt = ns_defaults()
%NS_DEFAULTS  Default option struct for NS_SOLVER.
%
%   opt = ns_defaults();  then override fields and call ns_solver(opt).
%
%   GRID / DOMAIN
%     nx, ny      number of CELLS                        [64 64]
%     Lx, Ly      domain size [m]                        [1e-3 1e-3]
%
%   FLUID
%     rho         density [kg/m^3]                       998
%     nu          kinematic viscosity [m^2/s]            1e-6
%     gx, gy      body acceleration [m/s^2]              0, 0
%     dpdx        applied mean pressure gradient [Pa/m]  0
%                 (used for periodic channel driving)
%
%   TIME
%     tEnd        physical end time [s]                  1e-3
%     cfl         safety factor on the timestep          0.5
%     dtMax       hard cap on dt [s]                     Inf
%     steady      stop early at steady state             false
%     steadyTol   tolerance on the residual              1e-4
%                 resid = max|u^{n+1}-u^n| / dt / Uref, a RATE, so the same
%                 tolerance means the same thing at any dt or grid size
%
%   NUMERICS
%     incremental   carry p^n in the predictor (2nd order)   true
%     implicitVisc  Crank-Nicolson on the viscous term       false
%     report        print diagnostics every N steps (0=off)  200
%
%   BOUNDARY CONDITIONS   opt.bc.<side> for side = left,right,bottom,top
%     .type = 'wall'      no-slip stationary wall
%           = 'moving'    no-slip wall moving tangentially (.u or .v)
%           = 'slip'      free slip (zero shear, no penetration)
%           = 'inflow'    prescribed normal velocity profile (.u or .v)
%           = 'outflow'   zero-gradient velocity, p = 0
%           = 'periodic'  paired with the opposite side
%
%   The default is a lid-driven cavity at 1 mm/s: three no-slip walls and a
%   moving top.  That is the Tutorial 7 validation case.
%
%   See also NS_SOLVER, CAVITY, CAVITY_VALIDATE.

opt.nx = 64;      opt.ny = 64;
opt.Lx = 1e-3;    opt.Ly = 1e-3;

opt.rho = 998;    opt.nu = 1e-6;
opt.gx  = 0;      opt.gy = 0;
opt.dpdx = 0;

opt.tEnd      = 1e-3;
opt.cfl       = 0.5;
opt.dtMax     = Inf;
opt.steady    = false;
opt.steadyTol = 1e-4;

opt.incremental  = true;
opt.implicitVisc = false;
opt.report       = 200;

% stokes: drop the convective term entirely.  Justified whenever Re << 1
% (Tutorial 1), and it also removes the cell-Reynolds instability that central
% differencing of convection suffers once Re_cell = |u|*dx/nu exceeds ~2.
% Recommended for any microfluidic geometry with corners.  Leave FALSE for the
% cavity benchmarks, which need inertia to be a meaningful test.
opt.stokes = false;

opt.bc.left   = struct('type','wall');
opt.bc.right  = struct('type','wall');
opt.bc.bottom = struct('type','wall');
opt.bc.top    = struct('type','moving','u',1e-3);

opt.u0 = [];      % optional initial fields (else zero)
opt.v0 = [];

% Spatially varying body force on the INTERIOR faces, used by the method of
% manufactured solutions (Tutorial 13).  Sizes: fx is (nx-1)-by-ny, fy is
% nx-by-(ny-1).  Leave empty for ordinary runs.
opt.fx = [];
opt.fy = [];

% Solid mask: nx-by-ny logical, true inside a solid body.  Velocity is forced
% to zero on every face touching a solid cell (direct-forcing immersed
% boundary).  See MAKE_MASK and the limitations in Tutorial 8.3 -- the boundary
% is staircased, so forces on the body are first-order accurate at best.
opt.solid = [];
end
