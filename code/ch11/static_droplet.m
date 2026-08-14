%STATIC_DROPLET  The standard two-phase benchmark: spurious currents.
%
%   Tutorial 11.  Run:  static_droplet
%
%   A droplet at rest, no gravity, no flow.  EXACTLY, nothing should happen:
%   surface tension is balanced by the pressure jump and u = 0 forever.
%
%   Numerically you get small vortices at the interface that never decay,
%   because the discrete sigma*kappa*delta*grad(phi) and the discrete grad(p)
%   are computed by different stencils and do not cancel.  These are SPURIOUS
%   (parasitic) currents, and how small they are is the single best measure of
%   a two-phase code's quality.
%
%   Checks:
%     1. pressure jump = sigma/R    (2-D: a cylinder, NOT 2*sigma/R)
%     2. spurious velocity small compared with the capillary scale sigma/mu
%     3. spurious velocity not GROWING
%     4. area (mass) conserved
%
%   See also NS2PHASE, LEVEL_SET, CURVATURE, YOUNG_LAPLACE.

if ~exist('QUIET','var'), QUIET = false; end

R  = 40e-6;
L  = 200e-6;
Ns = [32 64];

fprintf('\n  ===== Static droplet: spurious currents =====\n');
fprintf('  R = %g um in a %g um box. Exact answer: u = 0 everywhere,\n', R*1e6, L*1e6);
fprintf('  dp = sigma/R (2-D cylinder).\n');

res = struct([]);
for q = 1:numel(Ns)
    N = Ns(q);
    opt = ns2_defaults();
    opt.nx = N;  opt.ny = N;  opt.Lx = L;  opt.Ly = L;
    G = mac_grid(N,N,L,L);
    opt.phi0 = level_set('circle', G, L/2, L/2, R);
    opt.gx = 0;  opt.gy = 0;
    opt.Uref = 1e-4;
    opt.report = 0;

    % Run for a fixed number of capillary times so the comparison across grids
    % is like-for-like.
    tcap = sqrt((opt.rho1+opt.rho2)*R^3/(4*pi*opt.sigma));
    opt.tEnd = 60*tcap;

    fprintf('\n  --- N = %d (%.1f cells per radius) ---\n', N, R/(L/N));
    tic; S = ns2phase(opt); el = toc;

    sigma = opt.sigma;
    mu2   = opt.mu2;
    Ucap  = sigma/mu2;                       % capillary velocity scale

    % Pressure jump: inside mean minus outside mean, well away from the smear.
    inMask  = S.phi >  2*S.eps;
    outMask = S.phi < -2*S.eps;
    dp_num  = mean(S.p(inMask)) - mean(S.p(outMask));
    dp_ex   = sigma/R;

    umax = max(abs([S.u(:); S.v(:)]));
    % Growth check: compare the last quarter of the run with the middle.
    h = S.hist.umax;  m = numel(h);
    growth = mean(h(round(0.75*m):end)) / mean(h(round(0.4*m):round(0.6*m)));

    fprintf('  %d steps in %.1f s (dt limited by %s)\n', S.steps, el, S.dtLimit);
    fprintf('  pressure jump: numeric %8.2f Pa | exact sigma/R %8.2f Pa | err %5.1f%%\n', ...
            dp_num, dp_ex, 100*abs(dp_num-dp_ex)/dp_ex);
    fprintf('    (2*sigma/R would be %8.2f Pa -- the sphere formula, WRONG in 2-D)\n', ...
            2*sigma/R);
    fprintf('  spurious |u|max = %.3e m/s = %.2e * (sigma/mu)\n', umax, umax/Ucap);
    fprintf('  growth over the run (late/mid) = %.3f  (want <= ~1)\n', growth);
    fprintf('  area error = %+.3e\n', S.areaErr);

    res(q).N=N; res(q).umax=umax; res(q).Ucap=Ucap; res(q).dp=dp_num;
    res(q).dpex=dp_ex; res(q).area=S.areaErr; res(q).S=S; res(q).growth=growth;
end

fprintf('\n  Summary\n');
fprintf('     N   cells/R   |u|/(sigma/mu)   dp err     area err\n');
fprintf('   ---------------------------------------------------------\n');
for q = 1:numel(res)
    fprintf('   %4d   %6.1f     %12.2e   %6.1f%%   %+10.2e\n', ...
            res(q).N, R/(L/res(q).N), res(q).umax/res(q).Ucap, ...
            100*abs(res(q).dp-res(q).dpex)/res(q).dpex, res(q).area);
end

%% ------------------------------------------------------- THE VERDICT -----
% Do not leave this comparison to the reader.  A spurious-current number is
% meaningless in isolation; it only matters relative to the flow you intend to
% simulate.  So state the intended flow and judge against it.
Uphys = 1e-3;                             % a typical microchannel speed
uSpur = res(end).umax;
fprintf('\n  ===== VERDICT at N = %d =====\n', res(end).N);
fprintf('  spurious velocity      : %.2e m/s\n', uSpur);
fprintf('  intended physical flow : %.2e m/s\n', Uphys);
fprintf('  ratio                  : %.1f\n', uSpur/Uphys);
if uSpur < 0.1*Uphys
    fprintf('  PASS: parasitic flow is <10%% of the physical flow.\n');
else
    fprintf(['  *** FAIL for quantitative work at this grid. ***\n' ...
             '  The parasitic flow is LARGER than the flow you want to model,\n' ...
             '  so droplet trajectories and breakup times computed here would\n' ...
             '  be dominated by discretization error, not physics.\n' ...
             '\n  This is a real, reproducible limitation of this implementation\n' ...
             '  at these fluid properties (sigma/mu = %.2f m/s is a very fast\n' ...
             '  capillary scale). Options, in order of effectiveness:\n' ...
             '    1. Refine. The trend above shows it converging, but slowly.\n' ...
             '    2. Raise mu_c or lower sigma so sigma/mu falls -- often\n' ...
             '       physically legitimate (more viscous oil, more surfactant).\n' ...
             '    3. Use a sharper interface method (height functions for\n' ...
             '       curvature, or CLSVOF) -- the real fix, beyond this course.\n' ...
             '  Tutorial 12 results are therefore QUALITATIVE: trust the\n' ...
             '  regimes and trends, not the absolute droplet volumes.\n'], ...
             res(end).Ucap);
end
fprintf(['\n  Spurious currents never vanish and converge only slowly. That is\n' ...
         '  the honest state of the art, not a defect peculiar to this code --\n' ...
         '  but "everyone has this problem" is not a reason to skip measuring\n' ...
         '  YOUR value of it before quoting a result.\n\n']);

if ~QUIET
    S = res(end).S;  G = S.G;
    figure('Color','w','Position',[80 80 980 380]);
    subplot(1,3,1);
    contourf(G.xp*1e6,G.yp*1e6,S.p.',20,'LineStyle','none'); hold on
    contour(G.xp*1e6,G.yp*1e6,S.phi.',[0 0],'w','LineWidth',2);
    axis equal tight; colorbar; title('pressure [Pa]'); xlabel('x [\mum]');
    subplot(1,3,2);
    sp = sqrt(S.uc.^2+S.vc.^2);
    contourf(G.xp*1e6,G.yp*1e6,sp.',20,'LineStyle','none'); hold on
    contour(G.xp*1e6,G.yp*1e6,S.phi.',[0 0],'w','LineWidth',2);
    axis equal tight; colorbar; title('spurious |u| [m/s]');
    subplot(1,3,3);
    semilogy(S.hist.t*1e3, S.hist.umax,'LineWidth',1.5); grid on
    xlabel('t [ms]'); ylabel('max|u| [m/s]'); title('is it growing?');
end
