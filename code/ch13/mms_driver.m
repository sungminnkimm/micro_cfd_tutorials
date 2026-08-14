%MMS_DRIVER  Order-of-accuracy verification of ns_solver by manufactured solution.
%
%   Tutorial 13.  Run:  mms_driver
%
%   THE STRONGEST VERIFICATION TEST THERE IS.  Conservation and symmetry checks
%   can pass with a bug that happens to preserve the tested property.  An
%   order-of-accuracy test will not: a first-order bug ANYWHERE in a
%   second-order code drags the whole observed rate to 1, reliably.
%
%   THE MANUFACTURED SOLUTION
%   -------------------------
%   Derived from a streamfunction so it is divergence-free BY CONSTRUCTION:
%
%       psi = A sin^2(pi x) sin^2(pi y)
%       u   =  d(psi)/dy =  a sin^2(pi x) sin(2 pi y)
%       v   = -d(psi)/dx = -a sin(2 pi x) sin^2(pi y)      with a = A pi
%
%   Both components vanish on ALL FOUR walls, so the boundary condition is
%   plain no-slip everywhere -- no special treatment, and the test therefore
%   exercises the same BC code path as a real run.
%
%   With p_exact = 0, substituting into steady Navier-Stokes gives the required
%   body force
%
%       f = (u.grad)u - nu*lap(u)
%
%   which is added through opt.fx / opt.fy.  The solver must then reproduce u
%   exactly, and the error must fall at second order.
%
%   See also GCI, NS_SOLVER, POISSON_MMS.

if ~exist('QUIET','var'), QUIET = false; end

L  = 1;
nu = 0.05;             % large enough that the run reaches steady state quickly
rho = 1;
A  = 0.05;
a  = A*pi;  s = pi;

% --- exact fields and the derivatives needed for the forcing --------------
uex   = @(x,y)  a*sin(s*x).^2 .* sin(2*s*y);
vex   = @(x,y) -a*sin(2*s*x) .* sin(s*y).^2;
ux    = @(x,y)  a*s*sin(2*s*x).*sin(2*s*y);
uy    = @(x,y)  2*a*s*sin(s*x).^2 .* cos(2*s*y);
uxx   = @(x,y)  2*a*s^2*cos(2*s*x).*sin(2*s*y);
uyy   = @(x,y) -4*a*s^2*sin(s*x).^2 .* sin(2*s*y);
vx    = @(x,y) -2*a*s*cos(2*s*x).*sin(s*y).^2;
vy    = @(x,y) -a*s*sin(2*s*x).*sin(2*s*y);
vxx   = @(x,y)  4*a*s^2*sin(2*s*x).*sin(s*y).^2;
vyy   = @(x,y) -2*a*s^2*sin(2*s*x).*cos(2*s*y);

fprintf('\n  ===== MMS verification of ns_solver =====\n');
fprintf('  divergence-free by construction; check: max|div u_ex| on a fine grid\n');
nchk = 200; hchk = L/nchk;
xchk = (hchk/2:hchk:L-hchk/2)'; [Xk,Yk] = ndgrid(xchk,xchk);
fprintf('    max|du/dx + dv/dy| = %.3e   (analytic, should be ~0)\n', ...
        max(abs(ux(Xk,Yk)+vy(Xk,Yk)),[],'all'));

Ns  = [16 24 32 48 64];
errU = zeros(size(Ns));  errV = errU;  hh = errU;

for q = 1:numel(Ns)
    N = Ns(q);
    G = mac_grid(N,N,L,L);
    hh(q) = G.dx;

    opt = ns_defaults();
    opt.nx=N; opt.ny=N; opt.Lx=L; opt.Ly=L;
    opt.rho=rho; opt.nu=nu;
    opt.bc.left=struct('type','wall');  opt.bc.right=struct('type','wall');
    opt.bc.bottom=struct('type','wall');opt.bc.top=struct('type','wall');

    % Forcing on the INTERIOR faces, evaluated analytically at face centres.
    Xu = G.Xu(2:end-1,:);  Yu = G.Yu(2:end-1,:);
    Xv = G.Xv(:,2:end-1);  Yv = G.Yv(:,2:end-1);

    opt.fx = uex(Xu,Yu).*ux(Xu,Yu) + vex(Xu,Yu).*uy(Xu,Yu) ...
             - nu*(uxx(Xu,Yu) + uyy(Xu,Yu));
    opt.fy = uex(Xv,Yv).*vx(Xv,Yv) + vex(Xv,Yv).*vy(Xv,Yv) ...
             - nu*(vxx(Xv,Yv) + vyy(Xv,Yv));

    opt.Uref = max(abs(uex(Xk,Yk)),[],'all');
    opt.tEnd = 40;
    opt.steady = true;  opt.steadyTol = 1e-9;
    opt.report = 0;

    S = ns_solver(opt);

    eU = S.u - uex(G.Xu, G.Yu);
    eV = S.v - vex(G.Xv, G.Yv);
    errU(q) = norm(eU(:), inf);
    errV(q) = norm(eV(:), inf);

    fprintf('  N = %3d: %6d steps, max|div u| = %.2e, errU = %.3e, errV = %.3e\n', ...
            N, S.steps, S.divMax, errU(q), errV(q));
end

pU = log(errU(1:end-1)./errU(2:end)) ./ log(hh(1:end-1)./hh(2:end));
pV = log(errV(1:end-1)./errV(2:end)) ./ log(hh(1:end-1)./hh(2:end));

fprintf('\n     N      h        err(u)     order      err(v)     order\n');
fprintf('   ------------------------------------------------------------\n');
for q = 1:numel(Ns)
    fprintf('   %4d  %7.4f  %10.3e', Ns(q), hh(q), errU(q));
    if q>1, fprintf('  %7.2f', pU(q-1)); else, fprintf('      ---'); end
    fprintf('  %10.3e', errV(q));
    if q>1, fprintf('  %7.2f\n', pV(q-1)); else, fprintf('      ---\n'); end
end
fprintf('\n   observed order: u %.2f, v %.2f   (expected 2.00)\n', pU(end), pV(end));

fprintf(['\n   IF YOU MEASURE 1.0 instead of 2.0, the bug is almost always a\n' ...
         '   boundary condition -- see Exercise 13.2, which asks you to break\n' ...
         '   one deliberately and watch the rate fall. That is a much more\n' ...
         '   reliable bug detector than reading the code.\n\n']);

if ~QUIET
    figure('Color','w','Position',[80 80 620 460]);
    loglog(hh, errU,'o-','LineWidth',1.6); hold on
    loglog(hh, errV,'s-','LineWidth',1.6);
    loglog(hh, errU(1)*(hh/hh(1)).^2,'k--');
    grid on; xlabel('h'); ylabel('max-norm error');
    legend({'u','v','O(h^2)'},'Location','southeast');
    title('MMS: order of accuracy of the NS solver');
end
