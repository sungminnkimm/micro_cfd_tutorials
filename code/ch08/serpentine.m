%SERPENTINE  A snaking channel: the workhorse passive-micromixer geometry.
%
%   Tutorial 8.  Run:  serpentine
%
%   Produces the flow field consumed by Tutorial 9's mixing study, and saves it
%   to serpentine_flow.mat so you do not have to re-solve it every time.
%
%   WHAT 2-D CANNOT TELL YOU: at higher Re, bends generate DEAN VORTICES --
%   counter-rotating secondary flow in the channel CROSS-SECTION, which mixes
%   very effectively.  That is a 3-D effect and this 2-D model cannot see it.
%   The Dean number is De = Re*sqrt(h/(2R)); below De ~ 10 there is essentially
%   no secondary flow, so at microfluidic Re the omission is harmless.  Above
%   it, a 2-D serpentine will badly UNDERPREDICT mixing.
%
%   See also NS_SOLVER, MAKE_MASK, MIXING_STUDY, CHANNEL.

if ~exist('QUIET','var'), QUIET = false; end

%% ------------------------------------------------------------ geometry ---
Lx    = 2000e-6;      % total length
Ly    = 400e-6;       % bounding-box height
halfw = 50e-6;        % half channel width
amp   = 90e-6;        % amplitude of the snake
lam   = 500e-6;       % wavelength
Umean = 1e-3;

nx = 400;  ny = 80;
G  = mac_grid(nx, ny, Lx, Ly);

yCen  = @(x) Ly/2 + amp*sin(2*pi*x/lam);
solid = make_mask(G, {'below', @(x) yCen(x)-halfw}, ...
                     {'above', @(x) yCen(x)+halfw});

fprintf('\n  Serpentine: %g x %g um, channel half-width %g um\n', ...
        Lx*1e6, Ly*1e6, halfw*1e6);
fprintf('  grid %dx%d -> %.1f cells across the channel (want >= 20)\n', ...
        nx, ny, 2*halfw/G.dy);
fprintf('  open fraction %.2f\n', 1-mean(solid(:)));

%% ------------------------------------------------------------- solve -----
rho = 998;  mu = 1.0e-3;  nu = mu/rho;
opt = ns_defaults();
opt.nx = nx;  opt.ny = ny;  opt.Lx = Lx;  opt.Ly = Ly;
opt.rho = rho;  opt.nu = nu;
opt.solid = solid;

% Inlet: parabola across the OPEN part of the inlet cross-section only.
yq   = G.yu;
yLo  = yCen(0)-halfw;  yHi = yCen(0)+halfw;
uIn  = zeros(size(yq));
ins  = yq > yLo & yq < yHi;
eta  = (yq(ins)-yLo)/(2*halfw);
uIn(ins) = 6*Umean*eta.*(1-eta);

opt.bc.left   = struct('type','inflow','u',uIn);
opt.bc.right  = struct('type','outflow');
opt.bc.bottom = struct('type','wall');
opt.bc.top    = struct('type','wall');

% Re << 1 here, so inertia is negligible and the Stokes equations are the right
% model (Tutorial 1.2).  Dropping convection is also cheaper and removes one
% possible failure mode -- though NOT the one that actually bit this case; see
% Appendix B.4, where a wrong first diagnosis is worked through.
opt.stokes = true;
fprintf('  Re_cell = |u|dx/nu = %.4f (central differencing needs <~2: fine)\n', ...
        1.5*Umean*G.dx/nu);
fprintf('  -> solving the STOKES equations (opt.stokes = true), since Re << 1\n');

opt.tEnd   = 10*(2*halfw)^2/nu;
opt.steady = true;  opt.steadyTol = 1e-5;
opt.report = 0;

Re = Umean*2*halfw/nu;
fprintf('  Re = %.3g  |  De = Re*sqrt(h/2R) ~ %.3g (Dean flow needs De > 10)\n', ...
        Re, Re*sqrt(2*halfw/(2*lam)));

tic;  S = ns_solver(opt);  el = toc;
fprintf('\n  solved in %.1f s (%d steps)\n', el, S.steps);
fprintf('  mass balance %.3e | max|div u| %.3e | solid leak %.3e\n', ...
        S.massErr, S.divMax, S.solidLeak);

spd = sqrt(S.uc.^2 + S.vc.^2);
spd(solid) = NaN;
fprintf('  max speed %.4e m/s (%.2fx the mean -- bends accelerate the inside\n', ...
        max(spd(:),[],'omitnan'), max(spd(:),[],'omitnan')/Umean);
fprintf('  of the turn, which is where the extra shear comes from)\n');

%% ------------------------------------------------------------- save ------
outFile = fullfile(fileparts(mfilename('fullpath')),'serpentine_flow.mat');
save(outFile, 'S', 'G', 'solid', 'Umean', 'halfw', 'yCen', 'Lx', 'Ly');
fprintf('\n  saved flow field -> %s\n', outFile);
fprintf(['  Tutorial 9 loads this and transports a scalar on it. Because\n' ...
         '  Re << 1 the flow is STEADY, so solving it once and freezing it is\n' ...
         '  exact, not an approximation -- and it is what makes the mixing\n' ...
         '  study affordable.\n\n']);

if ~QUIET
    figure('Color','w','Position',[60 60 1100 340]);
    imagesc(G.xp*1e6, G.yp*1e6, spd.'*1e3);
    set(gca,'YDir','normal'); axis image; colorbar; hold on
    contour(G.xp*1e6, G.yp*1e6, double(solid.'), [0.5 0.5],'w','LineWidth',1.2);
    xlabel('x [\mum]'); ylabel('y [\mum]');
    title('Serpentine channel, |u| [mm/s]');
end
