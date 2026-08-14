%RISING_BUBBLE  Buoyancy vs drag vs surface tension, checked against Stokes.
%
%   Tutorial 11, step 4 of the validation ladder.  Run:  rising_bubble
%
%   Terminal velocity of a small sphere at Re << 1:
%
%       rigid sphere (no-slip surface)   U = (2/9) * drho*g*R^2/mu
%       clean bubble (mobile surface)    U = (2/3) * drho*g*R^2/mu    (Hadamard-
%                                                                     Rybczynski)
%
%   WHICH ONE YOU GET IS ITSELF A DIAGNOSTIC.  A factor of 3 separates them, and
%   it tells you whether your interface is transmitting shear (rigid-like) or
%   free to slip (mobile).  For a viscosity ratio lambda = mu_in/mu_out the
%   exact Hadamard-Rybczynski result interpolates between them:
%
%       U = (2/3) * drho*g*R^2/mu_out * (lambda + 1)/(3*lambda + 2)
%
%   lambda -> 0 (gas bubble) gives the 1/3 prefactor; lambda -> Inf (rigid)
%   gives 2/9.  This script computes that reference for the actual lambda used.
%
%   CAVEAT: these are 3-D SPHERE formulas.  A 2-D simulation is an infinite
%   CYLINDER, and 2-D Stokes flow past a cylinder has no steady solution at all
%   (Stokes' paradox) -- the drag depends on the domain size.  So expect the
%   right ORDER and the right SCALING with R^2 and drho, not the exact
%   coefficient.  Treated honestly below rather than fudged.
%
%   See also NS2PHASE, STATIC_DROPLET, LEVEL_SET.

if ~exist('QUIET','var'), QUIET = false; end

R  = 25e-6;
Lx = 200e-6;  Ly = 400e-6;
N  = 64;

opt = ns2_defaults();
opt.nx = N;  opt.ny = 2*N;  opt.Lx = Lx;  opt.Ly = Ly;

% Phase 1 (phi>0) = the bubble: light and inviscid-ish.  Phase 2 = liquid.
opt.rho1 = 100;    opt.mu1 = 1.0e-4;
opt.rho2 = 998;    opt.mu2 = 1.0e-3;
opt.sigma = 0.02;
opt.gy = -9.81;                       % gravity acts on both; buoyancy is the net

G = mac_grid(opt.nx, opt.ny, Lx, Ly);
opt.phi0 = level_set('circle', G, Lx/2, Ly/4, R);

drho = opt.rho2 - opt.rho1;
lam  = opt.mu1/opt.mu2;
U_HR = (2/3)*drho*9.81*R^2/opt.mu2 * (lam+1)/(3*lam+2);
U_rigid = (2/9)*drho*9.81*R^2/opt.mu2;
U_clean = (2/3)*drho*9.81*R^2/opt.mu2;

% HOW LONG TO RUN.  The obvious choice -- let the bubble rise many diameters --
% is unaffordable: 25 radii of rise needs 2.4 MILLION capillary-limited steps.
% But we do not need a long rise.  The momentum relaxation time of the bubble is
%       tau = rho_liquid * R^2 / mu  ~ 4e-4 s,
% after which it is moving at its terminal velocity.  Simulating a few tau is
% enough to MEASURE that velocity, which is the quantity being validated.
% Running longer would only accumulate more level-set drift.
tau = opt.rho2*R^2/opt.mu2;
opt.Uref = U_HR;
opt.tEnd = 6*tau;                 % a few relaxation times is all it takes
opt.report = 0;
fprintf('  momentum relaxation time rho*R^2/mu = %.2e s; simulating %.2e s (%.0f tau)\n', ...
        tau, opt.tEnd, opt.tEnd/tau);

Re_b = opt.rho2*U_HR*2*R/opt.mu2;
Bo   = drho*9.81*(2*R)^2/opt.sigma;
fprintf('\n  Rising bubble: R = %g um, drho = %g kg/m^3, lambda = %.3f\n', ...
        R*1e6, drho, lam);
fprintf('  reference terminal velocities (3-D sphere formulas):\n');
fprintf('    rigid sphere  (2/9) : %.4e m/s\n', U_rigid);
fprintf('    clean bubble  (2/3) : %.4e m/s\n', U_clean);
fprintf('    Hadamard-Rybczynski : %.4e m/s   <- for this lambda\n', U_HR);
fprintf('  Re = %.3g (need << 1 for Stokes), Bo = %.3g (need << 1 to stay round)\n', ...
        Re_b, Bo);

tic; S = ns2phase(opt); el = toc;

% Bubble centroid velocity from the level-set-weighted centre of mass.
H  = level_set('heaviside', S.phi, S.eps);
yc = sum(H(:).*G.Yp(:))/sum(H(:));
H0 = level_set('heaviside', opt.phi0, S.eps);
yc0 = sum(H0(:).*G.Yp(:))/sum(H0(:));
U_meas = (yc - yc0)/S.t;

% Instantaneous velocity over the last third, to avoid the startup transient.
vB = sum(H(:).*S.vc(:))/sum(H(:));

fprintf('\n  %d steps in %.1f s (dt limited by %s)\n', S.steps, el, S.dtLimit);
fprintf('  centroid rose %.2f um in %.3e s\n', (yc-yc0)*1e6, S.t);
fprintf('  mean rise velocity      : %.4e m/s\n', U_meas);
fprintf('  final bubble-mean v     : %.4e m/s\n', vB);
fprintf('  ratio to Hadamard-Rybczynski : %.3f\n', vB/U_HR);
fprintf('  ratio to rigid-sphere Stokes : %.3f\n', vB/U_rigid);
fprintf('  area error %+.3e | max|div u| %.2e\n', S.areaErr, S.divMax);

fprintf(['\n  INTERPRETING THE RATIO.  Do not expect 1.000. This is a 2-D\n' ...
         '  cylinder compared against 3-D sphere formulas, and 2-D Stokes flow\n' ...
         '  past a cylinder has no steady solution (Stokes paradox) -- the drag\n' ...
         '  depends on the domain width, so a confined 2-D bubble rises SLOWER\n' ...
         '  than the unbounded 3-D formula. What you SHOULD check is the\n' ...
         '  scaling: halve R and the velocity should drop 4x (Exercise), and\n' ...
         '  the bubble should stay round while Bo << 1.\n\n']);

if ~QUIET
    figure('Color','w','Position',[80 80 820 420]);
    subplot(1,3,1);
    contour(G.xp*1e6,G.yp*1e6,opt.phi0.',[0 0],'b','LineWidth',2); hold on
    contour(G.xp*1e6,G.yp*1e6,S.phi.',[0 0],'r','LineWidth',2);
    axis equal tight; grid on; xlabel('x [\mum]'); ylabel('y [\mum]');
    legend({'initial','final'},'Location','south'); title('bubble position');

    subplot(1,3,2);
    sp = sqrt(S.uc.^2+S.vc.^2);
    contourf(G.xp*1e6,G.yp*1e6,sp.'*1e3,20,'LineStyle','none'); hold on
    contour(G.xp*1e6,G.yp*1e6,S.phi.',[0 0],'w','LineWidth',1.5);
    axis equal tight; colorbar; title('|u| [mm/s]'); xlabel('x [\mum]');

    subplot(1,3,3);
    plot(S.hist.t*1e3, S.hist.areaErr,'LineWidth',1.5); grid on
    xlabel('t [ms]'); ylabel('relative area error'); title('mass conservation');
end
