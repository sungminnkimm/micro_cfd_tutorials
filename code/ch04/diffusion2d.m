%DIFFUSION2D  2-D unsteady diffusion, explicit vs implicit, with timing.
%
%   Tutorial 4, Exercise 4.6.  Run:  diffusion2d
%
%   The 2-D explicit limit is dt <= 1/(2D) * (1/dx^2 + 1/dy^2)^-1, which on a
%   uniform grid is dx^2/(4D) -- HALF the 1-D limit.  Forgetting to tighten dt
%   when moving from 1-D to 2-D is a rite of passage.
%
%   See also DIFFUSION1D, STABILITY_DEMO, POISSON2D.

if ~exist('QUIET','var'), QUIET = false; end

W = 100e-6;  H = 100e-6;  D = 5e-10;
tEnd = 0.15 * W^2/D;                     % ~3 s

grids = [32 64 128];
tExp = zeros(size(grids));  tImp = tExp;
nExp = zeros(size(grids));  nImp = nExp;

fprintf('\n  2-D diffusion to t = %.2f s.  Explicit dt ~ dx^2; implicit dt fixed.\n\n', tEnd);
fprintf('     N     dt_exp      steps |     dt_imp   steps |  t_exp   t_imp   ratio\n');
fprintf('  ---------------------------------------------------------------------\n');

for g = 1:numel(grids)
    n  = grids(g);
    dx = W/n;  dy = H/n;
    xc = (dx/2:dx:W-dx/2)';  yc = (dy/2:dy:H-dy/2)';
    [X,Y] = ndgrid(xc,yc);
    c0 = double((X-W/2).^2 + (Y-H/2).^2 < (W/6)^2);    % a disc of dye

    % ---------------- explicit ----------------
    dtE = 0.4 / (2*D) / (1/dx^2 + 1/dy^2);   % 80% of the 2-D limit
    nE  = ceil(tEnd/dtE);  dtE = tEnd/nE;
    rx = D*dtE/dx^2;  ry = D*dtE/dy^2;

    c = c0;  tic
    for k = 1:nE
        cp = padarray_neumann(c);                    % zero-flux ghost ring
        c  = c + rx*(cp(1:end-2,2:end-1) - 2*c + cp(3:end,2:end-1)) ...
               + ry*(cp(2:end-1,1:end-2) - 2*c + cp(2:end-1,3:end));
    end
    tExp(g) = toc;  nExp(g) = nE;  cE = c;

    % ---------------- implicit (Crank-Nicolson) ----------------
    % dt chosen for ACCURACY and held FIXED as the grid refines.
    dtI = tEnd/200;  nI = 200;
    ex = ones(n,1);
    Dxx = spdiags([ex -2*ex ex],[-1 0 1],n,n)/dx^2;  Dxx(1,1)=-1/dx^2; Dxx(n,n)=-1/dx^2;
    Dyy = spdiags([ex -2*ex ex],[-1 0 1],n,n)/dy^2;  Dyy(1,1)=-1/dy^2; Dyy(n,n)=-1/dy^2;
    Lap = kron(speye(n),Dxx) + kron(Dyy,speye(n));

    A  = speye(n*n) - dtI*D/2*Lap;
    B  = speye(n*n) + dtI*D/2*Lap;
    tic
    dA = decomposition(A,'lu');            % factorize ONCE (dt is constant)
    c  = c0(:);
    for k = 1:nI
        c = dA\(B*c);
    end
    tImp(g) = toc;  nImp(g) = nI;  cI = reshape(c,n,n);

    fprintf('  %4d  %9.2e  %9d | %9.2e  %6d | %6.2fs %6.2fs  %5.1fx\n', ...
            n, dtE, nE, dtI, nI, tExp(g), tImp(g), tExp(g)/tImp(g));

    if g == numel(grids)
        % Report mass RELATIVE to the initial value: absolute masses here are
        % ~1e-9 (areas in m^2) and print as 0.000000, which looks alarming and
        % means nothing.  Relative drift is the quantity with a right answer.
        m0 = sum(c0(:));  mE = sum(cE(:));  mI = sum(cI(:));
        fprintf('\n  Agreement explicit vs implicit on N=%d: max|diff| = %.3e\n', ...
                n, max(abs(cE(:)-cI(:))));
        fprintf('  Mass drift (relative to initial):  explicit %+.2e   implicit %+.2e\n', ...
                (mE-m0)/m0, (mI-m0)/m0);
        fprintf(['  Both are at round-off: zero-flux walls make total mass an\n' ...
                 '  exact invariant of both schemes, so anything above ~1e-12\n' ...
                 '  here would mean a leaking boundary condition.\n']);
        Xs=X; Ys=Y; cEs=cE; cIs=cI; c0s=c0; xs=xc; ys=yc;
    end
end

fprintf(['\n  Explicit step count grows as N^2 while each step also costs N^2,\n' ...
         '  so total explicit work ~ N^4.  Implicit holds the step count fixed.\n' ...
         '  That is the whole argument for implicit diffusion.\n\n']);

if ~QUIET
    figure('Color','w','Position',[100 100 940 320]);
    subplot(1,3,1); contourf(xs*1e6,ys*1e6,c0s.',20,'LineStyle','none');
    axis equal tight; colorbar; title('t = 0'); xlabel('x [\mum]');
    subplot(1,3,2); contourf(xs*1e6,ys*1e6,cEs.',20,'LineStyle','none');
    axis equal tight; colorbar; title(sprintf('explicit, t = %.1f s',tEnd));
    subplot(1,3,3); contourf(xs*1e6,ys*1e6,cIs.',20,'LineStyle','none');
    axis equal tight; colorbar; title('Crank-Nicolson');

    figure('Color','w','Position',[100 100 560 420]);
    loglog(grids, tExp,'o-','LineWidth',1.6); hold on
    loglog(grids, tImp,'s-','LineWidth',1.6);
    loglog(grids, tExp(1)*(grids/grids(1)).^4,'k--');
    grid on; xlabel('N'); ylabel('runtime [s]');
    legend({'explicit','Crank-Nicolson','O(N^4)'},'Location','northwest');
    title('cost of the explicit timestep restriction');
end
