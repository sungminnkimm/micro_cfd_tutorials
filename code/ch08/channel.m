%CHANNEL  Straight microchannel: the validation you run BEFORE adding geometry.
%
%   Tutorial 8.  Run:  channel
%
%   Every new geometry should start life as a straight channel you have
%   verified, and only then get its obstacle/bend/contraction added.  Debugging
%   a serpentine mixer from scratch is miserable; debugging one change to a
%   working channel is easy.
%
%   Checks performed:
%     * outlet profile vs the analytic parabola
%     * u_max/u_mean = 3/2
%     * global mass balance (should be ~1e-15)
%     * pressure drop vs  dp = 12 mu ubar L / h^2
%
%   See also NS_SOLVER, POISEUILLE, OBSTACLE, SERPENTINE.

if ~exist('QUIET','var'), QUIET = false; end

%% --------------------------------------------------------------- setup ---
h     = 100e-6;          % channel height [m]
Lx    = 800e-6;          % length [m]
Umean = 1e-3;            % target mean velocity [m/s]
rho   = 998;  mu = 1.0e-3;  nu = mu/rho;

nx = 160;  ny = 40;

opt = ns_defaults();
opt.nx = nx;  opt.ny = ny;
opt.Lx = Lx;  opt.Ly = h;
opt.rho = rho;  opt.nu = nu;

G  = mac_grid(nx, ny, Lx, h);
yc = G.yu;                                   % u-face heights = cell centres
uIn = 6*Umean * (yc/h) .* (1 - yc/h);        % parabola with mean = Umean

opt.bc.left   = struct('type','inflow','u',uIn);
opt.bc.right  = struct('type','outflow');
opt.bc.bottom = struct('type','wall');
opt.bc.top    = struct('type','wall');

% Time to steady state ~ viscous diffusion across the channel, h^2/nu.
opt.tEnd     = 5*h^2/nu;
opt.steady   = true;   opt.steadyTol = 1e-5;
opt.report   = 0;

Re = Umean*h/nu;
fprintf('\n  Straight channel: %g um x %g um, Umean = %g mm/s, Re = %.3g\n', ...
        Lx*1e6, h*1e6, Umean*1e3, Re);
fprintf('  grid %dx%d, entrance length 0.06*Re*h = %.3g um (negligible)\n', ...
        nx, ny, 0.06*Re*h*1e6);

tic;  S = ns_solver(opt);  el = toc;
fprintf('  solved in %.1f s (%d steps, t = %.3e s)\n\n', el, S.steps, S.t);

%% -------------------------------------------------------------- checks ---
uOut  = S.u(end,:).';
uExact = 6*Umean*(yc/h).*(1 - yc/h);
errP  = norm(uOut - uExact, inf)/max(uExact);

uMean_num = mean(uOut);
ratio     = max(uOut)/uMean_num;

% Pressure drop: compare cross-section averages near inlet and outlet, well
% clear of both ends so the numerical outlet condition does not contaminate it.
i1 = round(0.15*nx);  i2 = round(0.85*nx);
p1 = mean(S.p(i1,:));  p2 = mean(S.p(i2,:));
dp_num = p1 - p2;
dp_ex  = 12*mu*Umean*(G.xp(i2)-G.xp(i1))/h^2;

fprintf('  outlet profile vs parabola : %.3e  (relative, max-norm)\n', errP);
fprintf('  u_max/u_mean               : %.5f   (exact 1.5)\n', ratio);
fprintf('  global mass balance error  : %.3e   (should be ~1e-15)\n', S.massErr);
fprintf('  max|div u|                 : %.3e\n', S.divMax);
fprintf('  pressure drop  numeric %.5f Pa   analytic %.5f Pa  (%.2f%% err)\n', ...
        dp_num, dp_ex, 100*abs(dp_num-dp_ex)/dp_ex);
fprintf(['\n  The pressure-drop check is the sharpest of these: it exercises\n' ...
         '  the pressure solve, the boundary conditions and the viscous term\n' ...
         '  together. If dp is right to a fraction of a percent, the channel\n' ...
         '  solver works.\n']);

%% ------------------------------------------------- contraction variant ---
fprintf('\n  --- Contraction: %g um -> %g um ---\n', h*1e6, h/2*1e6);
opt2 = opt;
xStep = Lx/2;
% Solid blocks pinching the channel to half height over the downstream half.
opt2.solid = make_mask(G, {'rect', xStep, Lx, 0, h/4}, ...
                          {'rect', xStep, Lx, 3*h/4, h});
opt2.tEnd  = 5*h^2/nu;
opt2.steady = true;  opt2.steadyTol = 1e-5;
S2 = ns_solver(opt2);

% Mean velocity in the narrow section (only over the open cells).
iN = round(0.9*nx);
open = ~opt2.solid(iN,:);
uN = 0.5*(S2.u(iN,:) + S2.u(iN+1,:));
uMeanNarrow = sum(uN(open))/sum(open);

fprintf('  u_mean upstream  %.4e m/s\n', Umean);
fprintf('  u_mean in throat %.4e m/s   (mass conservation predicts %.4e)\n', ...
        uMeanNarrow, 2*Umean);
fprintf('  ratio %.3f  (exactly 2 if the throat is exactly half height)\n', ...
        uMeanNarrow/Umean);
fprintf('  mass balance %.3e | div %.3e | solid leak %.3e\n', ...
        S2.massErr, S2.divMax, S2.solidLeak);
fprintf(['\n  All three are at round-off. That is only true because the mask is\n' ...
         '  built into the PRESSURE OPERATOR (a solid face is a zero-flux face),\n' ...
         '  not merely zeroed on the velocities afterwards. Zeroing velocities\n' ...
         '  after the projection instead gives div ~ 1e2 here -- the flow simply\n' ...
         '  cannot be made divergence-free around a body the Poisson solve does\n' ...
         '  not know about. See poisson2d.m.\n\n' ...
         '  What masking still costs you is GEOMETRIC, not conservative: the\n' ...
         '  boundary is staircased, so forces on a curved body are first-order\n' ...
         '  accurate at best (Tutorial 8.3).\n\n']);

if ~QUIET
    figure('Color','w','Position',[80 80 980 560]);
    subplot(3,1,1);
    imagesc(G.xp*1e6, G.yp*1e6, sqrt(S.uc.^2+S.vc.^2).'*1e3);
    set(gca,'YDir','normal'); axis image; colorbar
    xlabel('x [\mum]'); ylabel('y [\mum]'); title('|u| [mm/s], straight channel');

    subplot(3,1,2);
    plot(uExact*1e3, yc*1e6,'-','LineWidth',2); hold on
    plot(uOut*1e3,   yc*1e6,'o','MarkerSize',4);
    grid on; xlabel('u [mm/s]'); ylabel('y [\mum]');
    legend({'analytic','outlet'},'Location','best'); title('outlet profile');

    subplot(3,1,3);
    sp = sqrt(S2.uc.^2+S2.vc.^2).'*1e3;  sp(opt2.solid.') = NaN;
    imagesc(G.xp*1e6, G.yp*1e6, sp);
    set(gca,'YDir','normal'); axis image; colorbar
    xlabel('x [\mum]'); ylabel('y [\mum]'); title('|u| [mm/s], contraction');
end
