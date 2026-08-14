%POISEUILLE  Pressure-driven flow between parallel plates on the MAC grid.
%
%   Tutorial 6.  Run:  poiseuille
%
%   Solves  mu d2u/dy2 = -G  with no-slip walls, on the staggered grid, and
%   checks against the analytic parabola
%
%       u(y) = G/(2 mu) * (h^2/4 - y^2),   u_max/u_mean = 3/2
%
%   WHY THIS IS NOT EXACT (and it is tempting to expect that it would be)
%   ---------------------------------------------------------------------
%   The INTERIOR stencil has zero truncation error on a quadratic, so you might
%   expect the parabola to be reproduced to round-off.  It is not, and the
%   reason is the cell-centred Dirichlet ghost.  With the wall at y = -h/2 and
%   cell 1 centred at y_1 = -h/2 + dy/2, the ghost sits at y_0 = -h/2 - dy/2 and
%   the TRUE parabola value there is
%
%       u(y_0) = A(-h*dy/2 - dy^2/4)     while    -u_1 = A(-h*dy/2 + dy^2/4)
%
%   so the linear ghost relation u_0 = -u_1 is wrong by A*dy^2/2.  That is an
%   O(dy^2) boundary error, which makes the scheme globally SECOND ORDER but
%   not exact.  Expect ~1e-4 relative error at ny = 64 and a clean slope of 2
%   under refinement -- verified below.
%
%   The lesson generalizes: "my scheme is second order" is a statement about
%   the interior AND the boundaries together, and the boundary is usually the
%   weaker of the two.  Always confirm the ORDER; only claim exactness when you
%   have derived that the boundary treatment is exact too.

%#ok<*NOPTS>
%
%   See also MAC_GRID, DUCT_CORRECTION, NS_SOLVER.

if ~exist('QUIET','var'), QUIET = false; end

%% ------------------------------------------------------------- physics ---
mu = 1.0e-3;                 % Pa.s, water
h  = 100e-6;                 % m, plate separation
G  = 1200;                   % Pa/m, -dp/dx  (gives ~1 mm/s, see Tutorial 1)

%% ---------------------------------------------------------- the solve ----
% 1-D in y is enough: the flow is fully developed and x-independent.  u lives
% at cell-centre height on the MAC grid, so the unknowns are u(y_j) with
% y_j = (j-1/2)dy, and BOTH walls need ghost cells (see MAC_GRID).
ny = 64;
dy = h/ny;
y  = (dy/2 : dy : h-dy/2)' - h/2;         % centred on the channel, in [-h/2,h/2]

% No-slip via ghost cells: u_ghost = 2*u_wall - u_adjacent = -u_adjacent.
% That is exactly the 'dirichlet0' operator from Tutorial 2.
[~, Lyy] = diff_matrix(ny, dy, 'dirichlet0');

u = (mu*Lyy) \ (-G*ones(ny,1));

%% ------------------------------------------------------------- checks ----
uex   = G/(2*mu) * (h^2/4 - y.^2);
umax  = G*h^2/(8*mu);
umean = G*h^2/(12*mu);

err       = norm(u - uex, inf);
umean_num = mean(u);                       % uniform cells -> plain mean is exact
umax_num  = max(u);
Q_num     = sum(u)*dy;                     % flow rate per unit depth
Q_cubic   = G*h^3/(12*mu);

% Wall shear stress: du/dy at the wall, using the ghost relation u_g = -u_1.
tau_num = mu * (u(1) - (-u(1))) / dy;      % = 2*mu*u(1)/dy
tau_ex  = G*h/2;

fprintf('\n  Poiseuille flow: h = %g um, G = %g Pa/m, mu = %g Pa.s, ny = %d\n', ...
        h*1e6, G, mu, ny);
fprintf('\n  max |u - u_exact|      : %.3e m/s\n', err);
fprintf('  relative to u_max      : %.3e  (boundary-limited, O(dy^2))\n', err/umax);
fprintf('\n  u_max   numeric %.6e   analytic %.6e\n', umax_num, umax);
fprintf('  u_mean  numeric %.6e   analytic %.6e\n', umean_num, umean);
fprintf('  u_max/u_mean = %.6f   (the 3/2 fingerprint for parallel plates)\n', ...
        umax_num/umean_num);
fprintf('\n  Q per unit depth  numeric %.6e   cubic law %.6e   (%.2e rel)\n', ...
        Q_num, Q_cubic, abs(Q_num-Q_cubic)/Q_cubic);
fprintf('  wall shear  numeric %.6f Pa   analytic %.6f Pa\n', tau_num, tau_ex);

%% ----------------------------------------- grid convergence (the check) --
% The honest version of "is my boundary treatment right?".  The interior is
% exact on a quadratic, so ANY error here is the boundary, and it must fall at
% second order.  A slope of 1 would mean the ghost formula is wrong.
fprintf('\n  Grid convergence (error is boundary-dominated):\n');
fprintf('     ny     max err     rel err    order\n');
fprintf('   ----------------------------------------\n');
nys = [16 32 64 128 256];
ec  = zeros(size(nys));
for k = 1:numel(nys)
    nyk = nys(k);  dyk = h/nyk;
    yk  = (dyk/2 : dyk : h-dyk/2)' - h/2;
    [~,Lk] = diff_matrix(nyk,dyk,'dirichlet0');
    uk  = (mu*Lk)\(-G*ones(nyk,1));
    ec(k) = norm(uk - G/(2*mu)*(h^2/4 - yk.^2), inf);
end
pc = log2(ec(1:end-1)./ec(2:end));
for k = 1:numel(nys)
    if k==1, fprintf('   %5d   %9.3e   %9.3e     ---\n', nys(k), ec(k), ec(k)/umax);
    else,    fprintf('   %5d   %9.3e   %9.3e    %5.2f\n', nys(k), ec(k), ec(k)/umax, pc(k-1));
    end
end
fprintf('   observed order %.2f (expected 2.00)\n', pc(end));

%% ---------------------------------------------- the cubic law, verified --
fprintf('\n  Cubic law check, Q ~ h^3 at fixed G:\n');
hs = [25 50 100 200]*1e-6;
Qs = zeros(size(hs));
for k = 1:numel(hs)
    nyk = 64; dyk = hs(k)/nyk;
    [~,Lk] = diff_matrix(nyk,dyk,'dirichlet0');
    uk = (mu*Lk)\(-G*ones(nyk,1));
    Qs(k) = sum(uk)*dyk;
end
slope = polyfit(log(hs), log(Qs), 1);
fprintf('    h [um] :'); fprintf(' %8.0f', hs*1e6); fprintf('\n');
fprintf('    Q      :'); fprintf(' %8.2e', Qs);    fprintf('\n');
fprintf('    d(logQ)/d(logh) = %.4f   (exactly 3 for the cubic law)\n', slope(1));

%% ------------------------------------------- how badly does 2-D lie? -----
fprintf('\n  The 2-D lie: correction factor for a RECTANGULAR duct\n');
fprintf('    w/h    correction   2-D overestimates flow rate by\n');
fprintf('   ------------------------------------------------------\n');
for ar = [1 2 5 10 20]
    f = duct_correction(ar);
    fprintf('   %5g     %7.4f          %6.0f%%\n', ar, f, 100*(1/f - 1));
end
fprintf(['\n  For a SQUARE channel, 2-D overpredicts flow rate by >130%%: the two\n' ...
         '  side walls you deleted were doing half the work. 2-D is a good model\n' ...
         '  for a shallow WIDE channel (w/h >~ 10) and a bad one near square.\n' ...
         '  It is much better for flow STRUCTURE than for flow MAGNITUDE.\n\n']);

if ~QUIET
    figure('Color','w','Position',[100 100 900 380]);
    subplot(1,2,1);
    plot(uex*1e3, y*1e6,'-','LineWidth',2); hold on
    plot(u*1e3,   y*1e6,'o','MarkerSize',4);
    grid on; xlabel('u [mm/s]'); ylabel('y [\mum]');
    legend({'analytic','MAC grid'},'Location','best');
    title(sprintf('Poiseuille, err = %.1e', err));

    subplot(1,2,2);
    ars = logspace(0,1.5,100);  fc = arrayfun(@duct_correction, ars);
    semilogx(ars, fc,'LineWidth',2); grid on; hold on
    yline(0.95,'k--','95%');
    xlabel('aspect ratio w/h'); ylabel('\bar{u} / \bar{u}_{2D}');
    title('how much of the 2-D answer survives in 3-D');
    ylim([0 1.05]);
end
