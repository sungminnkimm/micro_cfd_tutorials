%OBSTACLE  Cylindrical post in a microchannel: Stokes symmetry vs inertia.
%
%   Tutorial 8.  Run:  obstacle
%
%   THE TEST: at Re << 1 the streamlines around a symmetric body are FORE-AFT
%   SYMMETRIC -- the picture upstream is a mirror image of downstream.  That is
%   the visual signature of Stokes reversibility (Tutorial 1), and it is a
%   strong check: any visible asymmetry at Re = 0.01 means the convective term
%   has a bug or is being fed a wrong velocity.
%
%   At Re ~ 40 standing recirculation bubbles appear behind the post.  The
%   transition is around Re ~ 5.
%
%   See also NS_SOLVER, MAKE_MASK, CHANNEL.

if ~exist('QUIET','var'), QUIET = false; end

h   = 200e-6;         % channel height
Lx  = 1200e-6;        % length
rad = 40e-6;          % post radius
xc  = Lx/3;   yc0 = h/2;

nx = 240;  ny = 40;
G  = mac_grid(nx, ny, Lx, h);
fprintf('\n  Post radius %g um spans %.1f cells (want >= 20)\n', ...
        rad*1e6, 2*rad/G.dx);

ReList = [0.01 40];
out = cell(1,2);

for k = 1:2
    Re = ReList(k);
    rho = 998;  mu = 1.0e-3;  nu = mu/rho;
    Umean = Re*nu/h;

    opt = ns_defaults();
    opt.nx = nx;  opt.ny = ny;  opt.Lx = Lx;  opt.Ly = h;
    opt.rho = rho;  opt.nu = nu;
    yq  = G.yu;
    uIn = 6*Umean*(yq/h).*(1-yq/h);
    opt.bc.left   = struct('type','inflow','u',uIn);
    opt.bc.right  = struct('type','outflow');
    opt.bc.bottom = struct('type','wall');
    opt.bc.top    = struct('type','wall');
    opt.solid  = make_mask(G, {'circle', xc, yc0, rad});
    opt.tEnd   = 8*h^2/nu;
    opt.steady = true;  opt.steadyTol = 1e-5;
    opt.report = 0;

    fprintf('\n  --- Re = %g (Umean = %.3e m/s) ---\n', Re, Umean);
    tic;  S = ns_solver(opt);  el = toc;
    fprintf('  %d steps, %.1f s | mass err %.2e | div %.2e | solid leak %.2e\n', ...
            S.steps, el, S.massErr, S.divMax, S.solidLeak);

    % Fore-aft symmetry: compare centreline speed at equal distances up/down.
    jMid = round(ny/2);
    spd  = sqrt(S.uc.^2 + S.vc.^2);
    ds   = (2:6)*rad;
    asym = zeros(size(ds));
    for q = 1:numel(ds)
        iUp = find(G.xp >= xc-ds(q), 1);
        iDn = find(G.xp >= xc+ds(q), 1);
        a = spd(iUp,jMid);  b = spd(iDn,jMid);
        asym(q) = abs(a-b)/max(a,b);
    end
    fprintf('  fore-aft asymmetry at d = 2r..6r: ');
    fprintf('%.4f ', asym);  fprintf('\n');
    fprintf('  mean asymmetry = %.4f\n', mean(asym));

    out{k} = S;  out{k}.asym = asym;  out{k}.Umean = Umean;
end

fprintf(['\n  At Re = 0.01 the asymmetry should be a fraction of a percent --\n' ...
         '  residual staircasing and outlet-proximity effects, not physics.\n' ...
         '  At Re = 40 it is large and obvious: inertia has broken the\n' ...
         '  time-reversal symmetry of the Stokes equations.\n\n']);

if ~QUIET
    figure('Color','w','Position',[80 80 980 520]);
    for k = 1:2
        S = out{k};
        subplot(2,1,k);
        sp = sqrt(S.uc.^2+S.vc.^2).';
        sp(S.opt.solid.') = NaN;
        imagesc(G.xp*1e6, G.yp*1e6, sp/max(sp(:),[],'omitnan'));
        set(gca,'YDir','normal'); axis image; colorbar; hold on
        contour(G.xp*1e6, G.yp*1e6, double(S.opt.solid.'), [0.5 0.5],'w','LineWidth',1.5);
        xlabel('x [\mum]'); ylabel('y [\mum]');
        title(sprintf('Re = %g, normalized |u|  (mean asymmetry %.4f)', ...
                      ReList(k), mean(S.asym)));
    end
end
