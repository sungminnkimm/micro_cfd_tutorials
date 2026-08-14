%POISSON_DEMO  Steady heat conduction in a chip cross-section.
%
%   Tutorial 3, Exercises 3.2/3.3.  Run:  poisson_demo
%
%   200 um x 100 um silicon cross-section.  Hot bottom wall, cold top wall,
%   insulated sides, optional uniform Joule source.  Both cases have analytic
%   solutions, so this doubles as a validation.
%
%   See also POISSON2D, POISSON_SOLVE, POISSON_MMS.

if ~exist('QUIET','var'), QUIET = false; end

W  = 200e-6;    H  = 100e-6;      % m
kS = 148;                          % W/(m K), silicon
Tb = 350;       Tt = 300;          % K, bottom / top wall

nx = 96;  ny = 48;
dx = W/nx; dy = H/ny;
xc = (dx/2:dx:W-dx/2)';  yc = (dy/2:dy:H-dy/2)';

% Insulated sides -> Neumann in x;  fixed walls -> Dirichlet in y.
S = poisson2d(nx,ny,dx,dy,{'N','N','D','D'});

%% Case 1: no source.  Exact solution is a pure linear gradient in y.
T1 = poisson_solve(S, zeros(nx,ny), {0, 0, Tb, Tt});

Tex1 = Tb + (Tt-Tb)*(yc/H);                 % analytic
e1   = max(abs(T1 - repmat(Tex1.',nx,1)),[],'all');
fprintf('\nCase 1 (conduction, no source)\n');
fprintf('  max error vs analytic linear profile : %.3e K\n', e1);
fprintf('  x-variation (should be ~0)           : %.3e K\n', ...
        max(max(T1,[],1)-min(T1,[],1)));

% Heat flux through the top wall, q = -k dT/dy.  One-sided at the wall using
% the ghost relation T_ghost = 2*Tt - T(end).
dTdy_top = (Tt - T1(:,end)) / (dy/2);
q_num    = -kS * mean(dTdy_top);
q_ex     = -kS * (Tt-Tb)/H;
fprintf('  wall flux numeric %.4g W/m^2, analytic %.4g W/m^2 (%.2f%% err)\n', ...
        q_num, q_ex, 100*abs(q_num-q_ex)/abs(q_ex));

%% Case 2: uniform volumetric source (a Joule-heated resistor layer).
qv = 5e9;                          % W/m^3
% del^2 T = -qv/k
T2 = poisson_solve(S, -(qv/kS)*ones(nx,ny), {0, 0, Tb, Tt});

% Analytic: T = Tb + (Tt-Tb)*y/H + (qv/(2k))*y*(H-y)
Tex2 = Tb + (Tt-Tb)*(yc/H) + (qv/(2*kS))*yc.*(H-yc);
e2   = max(abs(T2 - repmat(Tex2.',nx,1)),[],'all');

% The meaningful metric is the PARABOLIC BULGE above the conduction baseline,
% not max(T)-Tb: with Tt < Tb the hottest point is still the bottom wall, so
% max(T)-Tb would report ~0 no matter how strong the source is.
bulge_num = max(T2 - repmat((Tb + (Tt-Tb)*(yc/H)).', nx, 1), [], 'all');
bulge_ex  = qv*H^2/(8*kS);

fprintf('\nCase 2 (with %.1e W/m^3 source)\n', qv);
fprintf('  max error vs analytic parabola : %.3e K\n', e2);
fprintf('  bulge above linear baseline    : %.4f K numeric, %.4f K analytic\n', ...
        bulge_num, bulge_ex);

% Same source, same geometry, but stagnant WATER instead of silicon.
kW = 0.6;
fprintf('  ...same source in stagnant water (k=%.1f): bulge would be %.2f K\n', ...
        kW, qv*H^2/(8*kW));
fprintf(['  Silicon conducts ~250x better, so the source barely registers.\n' ...
         '  This is why chips spread heat through silicon and why cooling a\n' ...
         '  microchannel needs FLOW, not conduction -- Tutorial 9 adds it.\n']);

if ~QUIET
    figure('Color','w','Position',[100 100 940 340]);
    subplot(1,2,1);
    contourf(xc*1e6, yc*1e6, T1.', 20, 'LineStyle','none');
    axis equal tight; colorbar; xlabel('x [\mum]'); ylabel('y [\mum]');
    title('conduction only');
    subplot(1,2,2);
    contourf(xc*1e6, yc*1e6, T2.', 20, 'LineStyle','none');
    axis equal tight; colorbar; xlabel('x [\mum]'); ylabel('y [\mum]');
    title(sprintf('with q = %.0e W/m^3', qv));

    figure('Color','w','Position',[100 100 520 400]);
    plot(T1(round(nx/2),:), yc*1e6,'o','MarkerSize',4); hold on
    plot(Tex1, yc*1e6,'-','LineWidth',1.5);
    plot(T2(round(nx/2),:), yc*1e6,'s','MarkerSize',4);
    plot(Tex2, yc*1e6,'-','LineWidth',1.5);
    grid on; xlabel('T [K]'); ylabel('y [\mum]');
    legend({'numeric (no src)','analytic','numeric (src)','analytic'}, ...
           'Location','best');
    title('centreline profiles vs analytic');
end
fprintf('\n');
