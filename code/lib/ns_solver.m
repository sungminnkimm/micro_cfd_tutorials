function S = ns_solver(opt)
%NS_SOLVER  Incompressible Navier-Stokes by the projection method (MAC grid).
%
%   S = ns_solver(opt)     opt from NS_DEFAULTS
%
%   ALGORITHM (Chorin 1968, incremental / van Kan variant)
%     1. predictor   u* = u^n + dt*(-conv + nu*lap + f - grad(p^n)/rho)
%     2. enforce the NORMAL component of u* at every boundary
%     3. Poisson     lap(phi) = div(u*)/dt
%     4. corrector   u^{n+1} = u* - dt*grad(phi)
%     5. pressure    p^{n+1} = p^n + rho*phi      (incremental only)
%
%   Step 2 is the one everybody skips.  The derivation of the phi boundary
%   condition (dphi/dn = 0 at a wall) ASSUMES u*_n equals u^{n+1}_n there.  Skip
%   it and you get a solver that produces plausible flows with a slow mass leak
%   -- the bug that takes a week to find.
%
%   STAGGERED LAYOUT (see MAC_GRID)
%     u : (nx+1) x ny     v : nx x (ny+1)     p : nx x ny
%
%   OUTPUT S
%     .u .v .p     staggered fields
%     .uc .vc      cell-centred velocities (for plotting/streamlines)
%     .G           grid struct from MAC_GRID
%     .divMax      final max|div u|   -- SHOULD BE ~1e-14
%     .t .steps    final time, step count
%     .hist        struct of time histories (t, ke, divMax)
%     .opt         the options actually used
%
%   THE DIAGNOSTIC THAT MATTERS: S.divMax.  If it is at round-off, your
%   projection is working and any remaining error is in the predictor.  That
%   splits the debugging problem cleanly in half.
%
%   See also NS_DEFAULTS, MAC_GRID, POISSON2D, CAVITY, CHANNEL.

if nargin < 1, opt = ns_defaults(); end

nx = opt.nx;  ny = opt.ny;
G  = mac_grid(nx, ny, opt.Lx, opt.Ly);
dx = G.dx;  dy = G.dy;  nu = opt.nu;  rho = opt.rho;

%% ------------------------------------------------ pressure Poisson setup --
% Wall/inflow -> Neumann (dphi/dn = 0).  Outflow -> Dirichlet (phi = 0).
bcT = { pbc(opt.bc.left), pbc(opt.bc.right), pbc(opt.bc.bottom), pbc(opt.bc.top) };
% The mask MUST go into the Poisson operator, not just onto the velocities.
% A solid face is a zero-flux face; if the operator does not know that, the
% projection cannot make the flow divergence-free around the body and you get
% div ~ 1e2 instead of 1e-14.
P   = poisson2d(nx, ny, dx, dy, bcT, opt.solid);
zeroBC = {zeros(ny,1), zeros(ny,1), zeros(nx,1), zeros(nx,1)};

%% ------------------------------------------------------------ solid mask --
% Expand the cell mask onto the faces: a face is blocked if EITHER neighbouring
% cell is solid.  Done once here rather than every step.
uMask = false(nx+1, ny);  vMask = false(nx, ny+1);
if ~isempty(opt.solid)
    if ~isequal(size(opt.solid),[nx ny])
        error('ns_solver:mask','opt.solid must be %dx%d.', nx, ny);
    end
    sld = logical(opt.solid);
    uMask(1:end-1,:) = uMask(1:end-1,:) | sld;
    uMask(2:end,:)   = uMask(2:end,:)   | sld;
    vMask(:,1:end-1) = vMask(:,1:end-1) | sld;
    vMask(:,2:end)   = vMask(:,2:end)   | sld;
end

%% --------------------------------------------------------- initial state --
u = zeros(nx+1, ny);   if ~isempty(opt.u0), u = opt.u0; end
v = zeros(nx, ny+1);   if ~isempty(opt.v0), v = opt.v0; end
p = zeros(nx, ny);

%% ------------------------------------------------------------- timestep ---
Uref = ref_velocity(opt);
dt_conv = dx / max(Uref, eps);
dt_visc = 0.5 / nu / (1/dx^2 + 1/dy^2);
if opt.implicitVisc
    dt = opt.cfl * dt_conv;
else
    dt = opt.cfl * min(dt_conv, dt_visc);
end
dt = min([dt, opt.dtMax, opt.tEnd]);
nsteps = max(1, ceil(opt.tEnd/dt));
dt = opt.tEnd/nsteps;

%% ---------------------------------------- implicit viscous factorization --
% Built ONCE, which is only valid because dt is held FIXED for the whole run.
if opt.implicitVisc
    Lu = build_lap_u(opt, nx, ny, dx, dy);
    Lv = build_lap_v(opt, nx, ny, dx, dy);
    dAu = decomposition(speye(size(Lu,1)) - dt*nu/2*Lu, 'lu');
    dAv = decomposition(speye(size(Lv,1)) - dt*nu/2*Lv, 'lu');
    Bu  = speye(size(Lu,1)) + dt*nu/2*Lu;
    Bv  = speye(size(Lv,1)) + dt*nu/2*Lv;
end

%% ---------------------------------------------------------- diagnostics ---
hist.t = zeros(nsteps+1,1); hist.ke = hist.t; hist.divMax = hist.t;
hist.resid = hist.t;
resid = Inf;

if opt.report
    fprintf('\n  ns_solver: %dx%d cells, dt = %.3e s, %d steps (t_end = %.3e s)\n', ...
            nx, ny, dt, nsteps, opt.tEnd);
    fprintf('    CFL limit %.3e s | viscous limit %.3e s | %s viscous\n', ...
            dt_conv, dt_visc, tern(opt.implicitVisc,'IMPLICIT','explicit'));
    fprintf('    Re = %.4g (based on Uref = %.3g m/s, L = %.3g m)\n\n', ...
            Uref*min(opt.Lx,opt.Ly)/nu, Uref, min(opt.Lx,opt.Ly));
    fprintf('      step        t [s]      max|u|     max|div u|        KE      resid\n');
    fprintf('    ---------------------------------------------------------------------\n');
end

%% =============================================================== TIME LOOP =
t = 0;
for n = 1:nsteps

    uPrev = u;

    % ---- 1. ghost-padded arrays carrying the velocity BCs ----------------
    [ue, ve] = apply_vel_bc(u, v, opt);

    % ---- 2. convection (divergence form) ---------------------------------
    % Interpolate to where the product lives, multiply THERE, difference back.
    % Never interpolate a product.
    uc  = 0.5*(u(1:end-1,:) + u(2:end,:));          % u^2 lives at centres
    vc  = 0.5*(v(:,1:end-1) + v(:,2:end));          % v^2 lives at centres
    ucor = 0.5*(ue(:,1:end-1) + ue(:,2:end));       % (nx+1)x(ny+1) corners
    vcor = 0.5*(ve(1:end-1,:) + ve(2:end,:));       % (nx+1)x(ny+1) corners
    uv   = ucor .* vcor;

    convU = (uc(2:end,:).^2 - uc(1:end-1,:).^2)/dx ...
          + (uv(2:end-1,2:end) - uv(2:end-1,1:end-1))/dy;

    convV = (vc(:,2:end).^2 - vc(:,1:end-1).^2)/dy ...
          + (uv(2:end,2:end-1) - uv(1:end-1,2:end-1))/dx;

    % STOKES OPTION: drop inertia entirely.  Tutorial 1 showed that the whole
    % left-hand side of Navier-Stokes is negligible at Re << 1, so this is the
    % physically correct model for a microchannel, and it is cheaper.
    %
    % It also removes one possible failure mode: central differencing of
    % convection is unconditionally unstable on its own (Tutorial 5.2) and is
    % stabilized only by viscosity, which needs the CELL Reynolds number
    %        Re_cell = |u| dx / nu  <~ 2.
    % Worth checking -- but check it with actual arithmetic. For a 100 um
    % channel at 1 mm/s on a 5 um grid, Re_cell = 1.5e-3 * 5e-6 / 1e-6 = 0.0075,
    % nowhere near the limit.  (An earlier draft of this comment claimed 7.5,
    % off by 1000, and sent the author chasing the wrong bug for an hour. The
    % real one was the pressure-outlet boundary value below.)
    if opt.stokes
        convU = 0*convU;  convV = 0*convV;
    end

    % ---- 3. viscous ------------------------------------------------------
    lapU = (ue(3:end,2:end-1) - 2*ue(2:end-1,2:end-1) + ue(1:end-2,2:end-1))/dx^2 ...
         + (ue(2:end-1,3:end) - 2*ue(2:end-1,2:end-1) + ue(2:end-1,1:end-2))/dy^2;

    lapV = (ve(3:end,2:end-1) - 2*ve(2:end-1,2:end-1) + ve(1:end-2,2:end-1))/dx^2 ...
         + (ve(2:end-1,3:end) - 2*ve(2:end-1,2:end-1) + ve(2:end-1,1:end-2))/dy^2;

    % ---- 4. pressure gradient (incremental form only) --------------------
    if opt.incremental
        gpU = (p(2:end,:) - p(1:end-1,:))/dx/rho;
        gpV = (p(:,2:end) - p(:,1:end-1))/dy/rho;
    else
        gpU = 0;  gpV = 0;
    end

    % ---- 5. predictor ----------------------------------------------------
    uI = u(2:end-1,:);   vI = v(:,2:end-1);        % interior unknowns

    % opt.fx/fy are spatially varying body-force arrays on the INTERIOR faces,
    % used by the method of manufactured solutions (Tutorial 13).  Empty for
    % ordinary runs.
    rhsU = -convU + opt.gx - gpU + opt.dpdx/rho;
    rhsV = -convV + opt.gy - gpV;
    if ~isempty(opt.fx), rhsU = rhsU + opt.fx; end
    if ~isempty(opt.fy), rhsV = rhsV + opt.fy; end

    if opt.implicitVisc
        % Split the Laplacian: the homogeneous part goes implicit, the
        % boundary contribution (what the explicit stencil sees minus what the
        % homogeneous operator gives) stays explicit.  Neat, and it keeps the
        % BC handling in ONE place -- apply_vel_bc.
        lapUhom = reshape(Lu*uI(:), size(uI));
        lapVhom = reshape(Lv*vI(:), size(vI));
        bcU = lapU - lapUhom;
        bcV = lapV - lapVhom;

        rU = uI(:) + dt*(rhsU(:) + nu*bcU(:));
        rV = vI(:) + dt*(rhsV(:) + nu*bcV(:));
        usI = reshape(dAu\(Bu*uI(:) + (rU - uI(:))), size(uI));
        vsI = reshape(dAv\(Bv*vI(:) + (rV - vI(:))), size(vI));
    else
        usI = uI + dt*(rhsU + nu*lapU);
        vsI = vI + dt*(rhsV + nu*lapV);
    end

    us = u;  us(2:end-1,:) = usI;
    vs = v;  vs(:,2:end-1) = vsI;

    % ---- 6. enforce NORMAL components of u* (see header) -----------------
    [us, vs] = apply_normal_bc(us, vs, opt);

    % ---- 6b. solid mask (direct-forcing immersed boundary) ---------------
    % Applied to u* BEFORE the Poisson solve so the divergence source knows
    % about the body, and again after the correction below.  See the honest
    % limitations in Tutorial 8.3: the boundary is staircased and the second
    % masking slightly breaks div = 0 right at the solid faces.  S.solidLeak
    % reports how much, so the error is measured rather than assumed.
    if ~isempty(opt.solid)
        us(uMask) = 0;  vs(vMask) = 0;
    end

    % ---- 7. pressure Poisson --------------------------------------------
    divS = (us(2:end,:) - us(1:end-1,:))/dx + (vs(:,2:end) - vs(:,1:end-1))/dy;

    % PRESSURE-OUTLET BOUNDARY VALUE FOR phi.  Subtle and important.
    % In the INCREMENTAL scheme the pressure accumulates, p^{n+1} = p^n + rho*phi.
    % Setting phi = 0 at a pressure outlet therefore pins the CORRECTION to
    % zero, not the PRESSURE -- so whatever p has drifted to at the outlet just
    % stays there and keeps growing. The straight channel hides this by
    % reaching steady state first; the serpentine runs longer and blows up at
    % the outlet face with p ~ 1e154.
    % What we actually want is p^{n+1} = 0 at the outlet, which needs
    %       phi_wall = -p_wall / rho.
    phiBC = zeroBC;
    if opt.incremental
        if bcT{1}=='D', phiBC{1} = -p(1,:).'  /rho; end
        if bcT{2}=='D', phiBC{2} = -p(end,:).'/rho; end
        if bcT{3}=='D', phiBC{3} = -p(:,1)    /rho; end
        if bcT{4}=='D', phiBC{4} = -p(:,end)  /rho; end
    end
    phi = poisson_solve(P, divS/dt, phiBC);

    % ---- 8. corrector ----------------------------------------------------
    u = us;  v = vs;
    u(2:end-1,:) = us(2:end-1,:) - dt*(phi(2:end,:) - phi(1:end-1,:))/dx;
    v(:,2:end-1) = vs(:,2:end-1) - dt*(phi(:,2:end) - phi(:,1:end-1))/dy;

    % BOUNDARY faces with a DIRICHLET pressure condition must be corrected too.
    % The Poisson operator assumed a flux through them (ghost = -phi_c, hence
    % the -2c on the diagonal); if the velocity there is instead extrapolated,
    % the projection and the operator disagree and the last cell column keeps a
    % finite divergence.  The face sits half a cell from the centre, so the
    % gradient is (phi_wall - phi_c)/(dx/2) with phi_wall = 0.
    % Use the SAME wall value the Poisson solve was given, or the operator and
    % the corrector disagree again (the bug from Tutorial 8.1, one level up).
    if bcT{1}=='D', u(1,:)   = us(1,:)   - dt*(phiBC{1}.' - phi(1,:))  *2/dx; end
    if bcT{2}=='D', u(end,:) = us(end,:) - dt*(phiBC{2}.' - phi(end,:))*2/dx; end
    if bcT{3}=='D', v(:,1)   = vs(:,1)   - dt*(phiBC{3}   - phi(:,1))  *2/dy; end
    if bcT{4}=='D', v(:,end) = vs(:,end) - dt*(phiBC{4}   - phi(:,end))*2/dy; end

    if ~isempty(opt.solid)
        solidLeak = max([max(abs(u(uMask)),[],'includenan'), ...
                         max(abs(v(vMask)),[],'includenan'), 0]);
        u(uMask) = 0;  v(vMask) = 0;
    else
        solidLeak = 0;
    end

    % ---- 9. pressure update ---------------------------------------------
    if opt.incremental
        p = p + rho*phi;
    else
        p = rho*phi;
    end

    t = t + dt;

    % ---- 10. diagnostics -------------------------------------------------
    divv = (u(2:end,:) - u(1:end-1,:))/dx + (v(:,2:end) - v(:,1:end-1))/dy;
    if ~isempty(opt.solid)
        divv(logical(opt.solid)) = 0;   % solid cells hold no fluid to conserve
    end
    ke   = 0.5*sum(u(:).^2)*dx*dy + 0.5*sum(v(:).^2)*dx*dy;

    % Steady-state residual: max rate of change of u, normalized by Uref.
    % Rate-based (not per-step), so the same tolerance means the same thing
    % regardless of dt or grid -- unlike a raw |KE^{n+1} - KE^n| test.
    resid = max(abs(u(:)-uPrev(:))) / dt / Uref;

    hist.t(n+1) = t;  hist.ke(n+1) = ke;  hist.divMax(n+1) = max(abs(divv(:)));
    hist.resid(n+1) = resid;

    if ~all(isfinite(u(:))) || ~all(isfinite(v(:)))
        error('ns_solver:diverged', ...
              ['Solution diverged at step %d (t = %.3e). Check the timestep ' ...
               '(dt = %.3e, viscous limit %.3e) and the BC setup.'], ...
              n, t, dt, dt_visc);
    end

    if opt.report && (mod(n,opt.report)==0 || n==1 || n==nsteps)
        fprintf('    %6d   %10.4e  %10.4e   %10.3e   %10.4e  %9.2e\n', ...
                n, t, max(abs(u(:))), max(abs(divv(:))), ke, resid);
    end

    if opt.steady && n > 10 && resid <= opt.steadyTol
        if opt.report
            fprintf('    steady state reached at step %d (t = %.4e s, resid = %.2e)\n', ...
                    n, t, resid);
        end
        break
    end
end
hist.t = hist.t(1:n+1);          hist.ke     = hist.ke(1:n+1);
hist.divMax = hist.divMax(1:n+1); hist.resid = hist.resid(1:n+1);
%% ===========================================================================

S.u = u;  S.v = v;  S.p = p;  S.G = G;
S.uc = 0.5*(u(1:end-1,:) + u(2:end,:));      % cell-centred, for plotting
S.vc = 0.5*(v(:,1:end-1) + v(:,2:end));
S.divMax = max(abs(divv(:)));
S.t = t;  S.steps = n;  S.dt = dt;
S.hist = hist;  S.opt = opt;
S.solidLeak = solidLeak;      % residual velocity on masked faces before re-zeroing
S.uMask = uMask;  S.vMask = vMask;

% Global mass balance -- the check from Tutorial 8.1.  Should be ~1e-15.
qIn  = sum(u(1,:))*dy;
qOut = sum(u(end,:))*dy;
if abs(qIn) > eps
    S.massErr = abs(qIn-qOut)/abs(qIn);
else
    S.massErr = abs(qIn-qOut);
end

if opt.report
    fprintf('\n    final max|div u| = %.3e', S.divMax);
    if S.divMax < 1e-10
        fprintf('   <-- at round-off, projection is working\n\n');
    else
        fprintf('   <-- TOO LARGE, something is wrong\n\n');
    end
end
end

%% ======================================================= helper functions ==

function t = pbc(b)
%PBC  Pressure BC type implied by a velocity BC type.
if strcmp(b.type,'outflow'), t = 'D'; else, t = 'N'; end
end

function s = tern(c,a,b)
if c, s = a; else, s = b; end
end

function U = ref_velocity(opt)
%REF_VELOCITY  Characteristic speed, for sizing the timestep.
U = 0;
sides = {'left','right','bottom','top'};
for k = 1:4
    b = opt.bc.(sides{k});
    if isfield(b,'u'), U = max(U, max(abs(b.u(:)))); end
    if isfield(b,'v'), U = max(U, max(abs(b.v(:)))); end
end
if ~isempty(opt.u0), U = max(U, max(abs(opt.u0(:)))); end
if ~isempty(opt.v0), U = max(U, max(abs(opt.v0(:)))); end
if opt.dpdx ~= 0
    h = min(opt.Lx, opt.Ly);
    U = max(U, abs(opt.dpdx)*h^2/(12*opt.rho*opt.nu));   % Poiseuille estimate
end
if U == 0, U = eps; end
end

function [ue, ve] = apply_vel_bc(u, v, opt)
%APPLY_VEL_BC  Build ghost-padded velocity arrays carrying the BCs.
%
%   ue is (nx+1)x(ny+2): u already sits ON the x-walls, so only y-ghosts are
%   needed.  ve is (nx+2)x(ny+1) for the mirror reason.  See MAC_GRID -- which
%   boundaries need ghosts SWAPS between u and v, and getting that backwards is
%   the classic staggered-grid bug.

% ---- u: set the x-boundary faces directly, then pad in y ----
u(1,:)   = wall_normal(opt.bc.left,  'u', u(2,:),     u(1,:));
u(end,:) = wall_normal(opt.bc.right, 'u', u(end-1,:), u(end,:));

ue = zeros(size(u,1), size(u,2)+2);
ue(:,2:end-1) = u;
ue(:,1)   = tang_ghost(opt.bc.bottom, 'u', u(:,1));
ue(:,end) = tang_ghost(opt.bc.top,    'u', u(:,end));

% ---- v: set the y-boundary faces directly, then pad in x ----
v(:,1)   = wall_normal(opt.bc.bottom, 'v', v(:,2),     v(:,1));
v(:,end) = wall_normal(opt.bc.top,    'v', v(:,end-1), v(:,end));

ve = zeros(size(v,1)+2, size(v,2));
ve(2:end-1,:) = v;
ve(1,:)   = tang_ghost(opt.bc.left,  'v', v(1,:));
ve(end,:) = tang_ghost(opt.bc.right, 'v', v(end,:));
end

function val = wall_normal(b, comp, inner, current)
%WALL_NORMAL  Value of the velocity component NORMAL to this wall (stored ON it).
switch b.type
    case {'wall','moving','slip'}
        val = zeros(size(inner));                 % no penetration
    case 'inflow'
        val = b.(comp) .* ones(size(inner));      % prescribed profile
        if isfield(b,comp) && numel(b.(comp))==numel(inner)
            val = reshape(b.(comp), size(inner));
        end
    case 'outflow'
        val = inner;                              % zero gradient
    case 'periodic'
        val = current;                            % handled by the caller
    otherwise
        error('ns_solver:bc','Unknown BC type ''%s''.', b.type);
end
end

function g = tang_ghost(b, comp, inner)
%TANG_GHOST  Ghost value for the velocity component TANGENTIAL to this wall.
%   No-slip:  u_ghost = 2*u_wall - u_inner   (wall value at the midpoint)
%   Slip:     u_ghost = u_inner              (zero shear)
switch b.type
    case 'wall'
        g = -inner;
    case 'moving'
        if isfield(b,comp), w = b.(comp); else, w = 0; end
        g = 2*w - inner;
    case {'slip','outflow'}
        g = inner;
    case 'inflow'
        g = -inner;                               % no-slip on the tangential part
    case 'periodic'
        g = inner;
    otherwise
        error('ns_solver:bc','Unknown BC type ''%s''.', b.type);
end
end

function [us, vs] = apply_normal_bc(us, vs, opt)
%APPLY_NORMAL_BC  Force u*_n = u^{n+1}_n at every boundary (step 2 of the header).
us(1,:)   = wall_normal(opt.bc.left,  'u', us(2,:),     us(1,:));
us(end,:) = wall_normal(opt.bc.right, 'u', us(end-1,:), us(end,:));
vs(:,1)   = wall_normal(opt.bc.bottom,'v', vs(:,2),     vs(:,1));
vs(:,end) = wall_normal(opt.bc.top,   'v', vs(:,end-1), vs(:,end));

% NOTE ON MASS BALANCE.  An earlier version rescaled the outflow profile here to
% force net flux to zero, on the reasoning that the Poisson equation needs the
% compatibility condition of Tutorial 3.  That reasoning applies to an
% ALL-NEUMANN problem.  With an outflow present, phi has a Dirichlet condition
% there, the operator is non-singular, no compatibility condition is required,
% and the projection sets the outlet flux itself -- which is exactly what a
% pressure outlet should do.  The zero-gradient assignment above is only an
% initial guess for u*; the pressure correction produces the final value, and
% global mass then balances to round-off as a CONSEQUENCE of div(u) = 0 rather
% than because it was imposed.
end

function L = build_lap_u(opt, nx, ny, dx, dy)
%BUILD_LAP_U  Homogeneous Laplacian on the (nx-1)-by-ny interior u unknowns.
%   x: u-faces sit ON the x-walls, so the interior unknowns are NODE-type with
%      known boundary values -> plain [1 -2 1] tridiagonal.
%   y: u sits at cell-centre height, walls half a cell outside -> CELL-centred,
%      so a no-slip wall gives -3 on the diagonal and a slip wall -1.
n  = nx-1;
ex = ones(n,1);
Dxx = spdiags([ex -2*ex ex],[-1 0 1],n,n)/dx^2;

ey = ones(ny,1);
Dyy = spdiags([ey -2*ey ey],[-1 0 1],ny,ny)/dy^2;
Dyy(1,1)   = Dyy(1,1)   + adj_cell(opt.bc.bottom)/dy^2;
Dyy(ny,ny) = Dyy(ny,ny) + adj_cell(opt.bc.top)/dy^2;

L = kron(speye(ny), Dxx) + kron(Dyy, speye(n));
end

function L = build_lap_v(opt, nx, ny, dx, dy)
%BUILD_LAP_V  Homogeneous Laplacian on the nx-by-(ny-1) interior v unknowns.
m  = ny-1;
ey = ones(m,1);
Dyy = spdiags([ey -2*ey ey],[-1 0 1],m,m)/dy^2;

ex = ones(nx,1);
Dxx = spdiags([ex -2*ex ex],[-1 0 1],nx,nx)/dx^2;
Dxx(1,1)   = Dxx(1,1)   + adj_cell(opt.bc.left)/dx^2;
Dxx(nx,nx) = Dxx(nx,nx) + adj_cell(opt.bc.right)/dx^2;

L = kron(speye(m), Dxx) + kron(Dyy, speye(nx));
end

function a = adj_cell(b)
%ADJ_CELL  Diagonal adjustment for a cell-centred tangential-velocity wall.
%   no-slip -> Dirichlet ghost (-1 adjustment, giving -3);
%   slip/outflow -> Neumann ghost (+1 adjustment, giving -1).
switch b.type
    case {'wall','moving','inflow'}, a = -1;
    case {'slip','outflow'},         a = +1;
    otherwise,                       a = -1;
end
end
