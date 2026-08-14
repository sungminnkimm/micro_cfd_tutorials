%FLOW_FOCUSING  Droplet generation by forcing two phases through an orifice.
%
%   Tutorial 12.  Run:  flow_focusing
%
%   Dispersed phase down the centre, continuous phase from both sides,
%   everything squeezed through a narrow orifice.  Compared with a T-junction
%   this gives SMALLER droplets (below the channel width), higher rates, and
%   finer control, because the ORIFICE width -- not the channel width -- sets
%   the scale:
%
%       d/w_orifice ~ (Qd/Qc)^(1/3) * Ca^(-1/3)     (dripping regime)
%
%   Same caveats as T_JUNCTION: 2-D gets trends right and absolute volumes
%   wrong, and the neck must be resolved or breakup becomes a numerical
%   artefact rather than a physical one.
%
%   See also T_JUNCTION, REGIME_MAP, NS2PHASE.

if ~exist('QUIET','var'), QUIET = false; end
if ~exist('RATIO','var'), RATIO = 0.5; end

%% ------------------------------------------------------------ geometry ---
%  Sized for affordability, exactly as in T_JUNCTION -- see the cautionary note
%  in that file about the 5.3-million-step version that came first.  Run
%  REGIME_MAP before changing any of these numbers.
Lx = 300e-6;  Ly = 150e-6;
wCh  = 90e-6;                  % upstream channel height (open region)
wOr  = 30e-6;                  % orifice height
xOr  = 130e-6;  LOr = 45e-6;   % orifice start and length

nx = 100;  ny = 50;
G  = mac_grid(nx, ny, Lx, Ly);

yMid = Ly/2;
% Channel walls: a straight duct of height wCh, pinched to wOr over the orifice.
halfH = @(x) (wCh/2)*ones(size(x)) - ((wCh-wOr)/2) * ...
             double(x >= xOr & x <= xOr+LOr);
solid = make_mask(G, {'above', @(x) yMid + halfH(x)}, ...
                     {'below', @(x) yMid - halfH(x)});

fprintf('\n  Flow focusing: orifice %g um at %.1f cells across\n', ...
        wOr*1e6, wOr/G.dy);

%% -------------------------------------------------------------- fluids ---
mu_c = 2.0e-3; rho_c = 900;
mu_d = 1.0e-3; rho_d = 998;
sigma = 5e-3;

Uc = 12e-3;
Ca = mu_c*Uc/sigma;
fprintf('  Ca = %.4f, Qd/Qc = %.2f\n', Ca, RATIO);
fprintf('  dripping-regime estimate d/w_or ~ (Qd/Qc)^(1/3)*Ca^(-1/3) = %.2f\n', ...
        (RATIO)^(1/3)*Ca^(-1/3));

opt = ns2_defaults();
opt.nx = nx;  opt.ny = ny;  opt.Lx = Lx;  opt.Ly = Ly;
opt.rho1 = rho_d;  opt.mu1 = mu_d;
opt.rho2 = rho_c;  opt.mu2 = mu_c;
opt.sigma = sigma;
opt = set_solid(opt, solid);

% Open boundaries: the dispersed volume is SUPPOSED to grow as fluid arrives,
% so the global area correction must be OFF here.
opt.conserveMass = false;

% Inlet: dispersed phase occupies the central third, continuous the rest.
yu   = G.yu;
open = abs(yu - yMid) < wCh/2;
wD   = wCh/3;
uIn  = zeros(size(yu));
prof = zeros(size(yu));
eta  = (yu(open) - (yMid-wCh/2))/wCh;
prof(open) = 6*Uc*eta.*(1-eta);
uIn = prof;
opt.bc.left   = struct('type','inflow','u',uIn);
opt.bc.right  = struct('type','outflow');
opt.bc.top    = struct('type','wall');
opt.bc.bottom = struct('type','wall');

% Initial interface: a tongue of dispersed phase along the centreline inlet.
phi0 = min(wD/2 - abs(G.Yp - yMid), 60e-6 - G.Xp);
phi0(solid) = -abs(phi0(solid));
opt.phi0 = phi0;

opt.Uref = Uc*3;
opt.tEnd = 0.6*Lx/Uc;
opt.report = 5000;

tic; S = ns2phase(opt); el = toc;
fprintf('  wall time %.1f s, %d steps, dt limited by %s\n', el, S.steps, S.dtLimit);

%% ------------------------------------------------------------ measure ----
jMid = round(yMid/G.dy);
line = S.phi(:,jMid) > 0;
downstream = G.xp > xOr + LOr;
seg = line & downstream;
[st, ln] = runs_true(seg, G.xp);

fprintf('\n  droplets downstream of the orifice: %d\n', numel(ln));
if ~isempty(ln)
    fprintf('    lengths [um]: '); fprintf('%.1f ', ln*1e6); fprintf('\n');
    fprintf('    mean d/w_or = %.2f\n', mean(ln)/wOr);
else
    fprintf('    None yet -- longer run, finer neck, or different Ca needed.\n');
    fprintf('    Check the diagnostics before treating this as physics.\n');
end
fprintf('\n  max|div u| %.2e | dt %.2e s\n\n', S.divMax, S.dt);

if ~QUIET
    figure('Color','w','Position',[60 60 1050 400]);
    ph = S.phi.'; ph(solid.') = NaN;
    imagesc(G.xp*1e6, G.yp*1e6, ph); set(gca,'YDir','normal'); hold on
    contour(G.xp*1e6, G.yp*1e6, S.phi.', [0 0],'k','LineWidth',2);
    contour(G.xp*1e6, G.yp*1e6, double(solid.'), [0.5 0.5],'w','LineWidth',1.5);
    axis image; colorbar; xlabel('x [\mum]'); ylabel('y [\mum]');
    title(sprintf('Flow focusing, Ca = %.4f', Ca));
end

function [starts, lens] = runs_true(mask, x)
mask = mask(:).';
d = diff([false mask false]);
i0 = find(d==1);  i1 = find(d==-1)-1;
starts = x(i0);
lens = x(min(i1,numel(x))) - x(i0) + (x(2)-x(1));
end
