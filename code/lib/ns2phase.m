function S = ns2phase(opt)
%NS2PHASE  Two-phase incompressible flow: projection + level set + CSF.
%
%   S = ns2phase(opt)    opt from NS2_DEFAULTS
%
%   Extends NS_SOLVER with:
%     * a level-set interface phi (phi > 0 = fluid 1 = the DISPERSED phase)
%     * variable density and viscosity, smeared over ~1.5 cells
%     * continuum surface force  f = sigma * kappa * delta(phi) * grad(phi)
%
%   BALANCED-FORCE DISCRETIZATION
%   -----------------------------
%   Spurious currents arise when the discrete surface-tension force and the
%   discrete pressure gradient do not cancel, because they are computed by
%   different stencils in different places.  The cure is to evaluate
%   sigma*kappa*delta*grad(phi) at the SAME FACES where grad(p) lives, using the
%   SAME difference operator:
%
%       f_x at u-face i  =  sigma * kappa_face * (phi(i) - phi(i-1))/dx
%
%   with kappa_face the average of the two adjacent cell curvatures.  Then a
%   static droplet's pressure jump can be represented EXACTLY by the discrete
%   pressure field, and the residual velocity drops by orders of magnitude
%   compared to a naive cell-centred force.  STATIC_DROPLET measures it.
%
%   VARIABLE DENSITY.  The pressure Poisson equation strictly becomes
%   div( grad(phi)/rho ) = div(u*)/dt.  This implementation uses the CONSTANT-
%   COEFFICIENT approximation with rho = min(rho1,rho2), which is standard for
%   modest density ratios and keeps the Tutorial 3 factorization reusable.  It
%   is accurate for water/oil (ratio ~1) and degrades for water/air (1000:1).
%   Start with water/oil.
%
%   TIMESTEP is the min of CFL, viscous, and CAPILLARY limits.  The solver
%   prints which one is binding -- more useful than memorizing the crossover.
%
%   See also NS2_DEFAULTS, LEVEL_SET, CURVATURE, STATIC_DROPLET, RISING_BUBBLE.

if nargin < 1, opt = ns2_defaults(); end

nx = opt.nx;  ny = opt.ny;
G  = mac_grid(nx, ny, opt.Lx, opt.Ly);
dx = G.dx;  dy = G.dy;
sig = opt.sigma;
ep  = opt.epsFac*max(dx,dy);

rho1 = opt.rho1;  rho2 = opt.rho2;
mu1  = opt.mu1;   mu2  = opt.mu2;
rhoP = min(rho1, rho2);            % constant-coefficient pressure density

%% ---------------------------------------------------- pressure operator --
bcT = { pbc(opt.bc.left), pbc(opt.bc.right), pbc(opt.bc.bottom), pbc(opt.bc.top) };
P   = poisson2d(nx, ny, dx, dy, bcT, opt.solid);
zbc = {zeros(ny,1), zeros(ny,1), zeros(nx,1), zeros(nx,1)};

%% --------------------------------------------------------- initial data --
phi = opt.phi0;
if isempty(phi), error('ns2phase:phi0','opt.phi0 is required.'); end
u = zeros(nx+1, ny);  if ~isempty(opt.u0), u = opt.u0; end
v = zeros(nx, ny+1);  if ~isempty(opt.v0), v = opt.v0; end
p = zeros(nx, ny);

%% ------------------------------------------------------------ timestep ---
Umax0 = max([max(abs(u(:))), max(abs(v(:))), opt.Uref, eps]);
nuMax = max(mu1/rho1, mu2/rho2);

dt_cfl  = dx/Umax0;
dt_visc = 0.5/nuMax/(1/dx^2 + 1/dy^2);
if sig > 0
    dt_cap = sqrt((rho1+rho2)*min(dx,dy)^3/(4*pi*sig));
else
    dt_cap = Inf;
end
[dtRaw, which] = min([dt_cfl, dt_visc, dt_cap]);
names = {'CFL','viscous','CAPILLARY'};
dt = opt.cfl*dtRaw;
dt = min([dt, opt.dtMax, opt.tEnd]);
nsteps = max(1, ceil(opt.tEnd/dt));
dt = opt.tEnd/nsteps;

A0 = level_set('area', phi, dx, dy, ep);

% GUARD.  A capillary- or viscous-limited two-phase run can silently demand
% millions of steps -- the first draft of T_JUNCTION asked for 5.3 million,
% about 90 minutes, and gave no warning until it failed to finish.  Refuse to
% start rather than discover that later.  Raise opt.maxSteps deliberately if
% you really do want a long run.
if nsteps > opt.maxSteps
    error('ns2phase:tooManySteps', ...
        ['This run needs %s steps (dt = %.2e s, limited by %s; t_end = %.3e s).\n' ...
         'That is above opt.maxSteps = %s. Options:\n' ...
         '  * coarsen the grid (dt_cap ~ dx^1.5, dt_visc ~ dx^2)\n' ...
         '  * lower mu of the continuous phase (raises dt_visc AND lowers Ca)\n' ...
         '  * raise the velocity to shorten the physical time simulated\n' ...
         '  * shorten opt.tEnd\n' ...
         'Run REGIME_MAP first next time -- it reports this in one second.'], ...
        num2str(nsteps), dt, names{which}, opt.tEnd, num2str(opt.maxSteps));
end

if opt.report
    fprintf('\n  ns2phase: %dx%d, sigma = %g N/m, eps = %.2f dx\n', ...
            nx, ny, sig, ep/dx);
    fprintf('    dt limits: CFL %.3e | viscous %.3e | capillary %.3e\n', ...
            dt_cfl, dt_visc, dt_cap);
    fprintf('    BINDING: %s   ->  dt = %.3e s, %d steps\n', names{which}, dt, nsteps);
    fprintf('    initial enclosed area = %.6e m^2\n\n', A0);
    fprintf('      step        t [s]      max|u|    max|div u|   area err   |grad phi|\n');
    fprintf('    -------------------------------------------------------------------\n');
end

hist.t = zeros(nsteps,1);  hist.umax = hist.t;  hist.areaErr = hist.t;
hist.divMax = hist.t;      hist.gmag = hist.t;

%% ================================================================= LOOP ==
t = 0;
for n = 1:nsteps

    %% ---- 1. properties from phi -------------------------------------
    H   = level_set('heaviside', phi, ep);
    rho = rho2 + (rho1-rho2)*H;
    mu  = mu2  + (mu1 -mu2 )*H;

    %% ---- 2. curvature (cell-centred) --------------------------------
    [kap, ~, ~, gmag] = curvature(phi, dx, dy);
    kap = max(min(kap, 1/dx), -1/dx);        % clip: |k| > 1/dx is unresolved

    dlt = level_set('delta', phi, ep);

    %% ---- 3. ghost-padded velocities ---------------------------------
    [ue, ve] = pad_vel(u, v, opt);

    %% ---- 4. convection ----------------------------------------------
    uc = 0.5*(u(1:end-1,:) + u(2:end,:));
    vc = 0.5*(v(:,1:end-1) + v(:,2:end));
    ucor = 0.5*(ue(:,1:end-1) + ue(:,2:end));
    vcor = 0.5*(ve(1:end-1,:) + ve(2:end,:));
    uv   = ucor.*vcor;

    convU = (uc(2:end,:).^2 - uc(1:end-1,:).^2)/dx ...
          + (uv(2:end-1,2:end) - uv(2:end-1,1:end-1))/dy;
    convV = (vc(:,2:end).^2 - vc(:,1:end-1).^2)/dy ...
          + (uv(2:end,2:end-1) - uv(1:end-1,2:end-1))/dx;

    %% ---- 5. viscous (variable mu, face-averaged) ---------------------
    lapU = (ue(3:end,2:end-1) - 2*ue(2:end-1,2:end-1) + ue(1:end-2,2:end-1))/dx^2 ...
         + (ue(2:end-1,3:end) - 2*ue(2:end-1,2:end-1) + ue(2:end-1,1:end-2))/dy^2;
    lapV = (ve(3:end,2:end-1) - 2*ve(2:end-1,2:end-1) + ve(1:end-2,2:end-1))/dx^2 ...
         + (ve(2:end-1,3:end) - 2*ve(2:end-1,2:end-1) + ve(2:end-1,1:end-2))/dy^2;

    muU = 0.5*(mu(1:end-1,:) + mu(2:end,:));        % (nx-1) x ny  at u-faces
    muV = 0.5*(mu(:,1:end-1) + mu(:,2:end));        % nx x (ny-1)  at v-faces
    rhoU = 0.5*(rho(1:end-1,:) + rho(2:end,:));
    rhoV = 0.5*(rho(:,1:end-1) + rho(:,2:end));

    %% ---- 6. BALANCED-FORCE surface tension ---------------------------
    % Evaluate at the SAME faces, with the SAME operator, as grad(p).
    kapU = 0.5*(kap(1:end-1,:) + kap(2:end,:));
    kapV = 0.5*(kap(:,1:end-1) + kap(:,2:end));
    dltU = 0.5*(dlt(1:end-1,:) + dlt(2:end,:));
    dltV = 0.5*(dlt(:,1:end-1) + dlt(:,2:end));

    fstU = sig * kapU .* dltU .* (phi(2:end,:) - phi(1:end-1,:))/dx;
    fstV = sig * kapV .* dltV .* (phi(:,2:end) - phi(:,1:end-1))/dy;

    %% ---- 7. predictor ------------------------------------------------
    gpU = (p(2:end,:) - p(1:end-1,:))/dx;
    gpV = (p(:,2:end) - p(:,1:end-1))/dy;

    usI = u(2:end-1,:) + dt*( -convU + (muU.*lapU + fstU - gpU)./rhoU + opt.gx );
    vsI = v(:,2:end-1) + dt*( -convV + (muV.*lapV + fstV - gpV)./rhoV + opt.gy );

    us = u;  us(2:end-1,:) = usI;
    vs = v;  vs(:,2:end-1) = vsI;
    [us, vs] = set_normal(us, vs, opt);
    if ~isempty(opt.solid)
        us(opt.uMask) = 0;  vs(opt.vMask) = 0;
    end

    %% ---- 8. projection -----------------------------------------------
    divS = (us(2:end,:)-us(1:end-1,:))/dx + (vs(:,2:end)-vs(:,1:end-1))/dy;

    % Pressure-outlet value for q.  This solver is incremental (p = p + q), so
    % q = 0 at an outlet pins the CORRECTION rather than the PRESSURE and lets
    % p drift without bound there.  See the same note in NS_SOLVER.  q already
    % carries pressure units, so the wall value is simply -p_wall.
    qbc = zbc;
    if bcT{1}=='D', qbc{1} = -p(1,:).';   end
    if bcT{2}=='D', qbc{2} = -p(end,:).'; end
    if bcT{3}=='D', qbc{3} = -p(:,1);     end
    if bcT{4}=='D', qbc{4} = -p(:,end);   end
    q = poisson_solve(P, rhoP*divS/dt, qbc);

    u = us;  v = vs;
    u(2:end-1,:) = us(2:end-1,:) - dt*(q(2:end,:)-q(1:end-1,:))/dx/rhoP;
    v(:,2:end-1) = vs(:,2:end-1) - dt*(q(:,2:end)-q(:,1:end-1))/dy/rhoP;
    if bcT{1}=='D', u(1,:)  = us(1,:)  - dt*(qbc{1}.'-q(1,:))  *2/dx/rhoP; end
    if bcT{2}=='D', u(end,:)= us(end,:)- dt*(qbc{2}.'-q(end,:))*2/dx/rhoP; end
    if bcT{3}=='D', v(:,1)  = vs(:,1)  - dt*(qbc{3}  -q(:,1))  *2/dy/rhoP; end
    if bcT{4}=='D', v(:,end)= vs(:,end)- dt*(qbc{4}  -q(:,end))*2/dy/rhoP; end
    if ~isempty(opt.solid), u(opt.uMask) = 0;  v(opt.vMask) = 0; end

    p = p + q;

    %% ---- 9. advect and reinitialize the interface ---------------------
    phi = level_set('advect', phi, u, v, dt, dx, dy);
    if mod(n, opt.reinitEvery) == 0
        phi = level_set('reinit', phi, dx, dy, opt.reinitIters);
        % Reinitialization drifts the zero contour by ~0.1% of the area per
        % call.  Over the hundreds of calls a droplet run needs, that alone
        % dissolves the droplet entirely.  Restore the target area unless the
        % problem has open boundaries, where the enclosed area is SUPPOSED to
        % change as fluid enters and leaves.
        if opt.conserveMass
            phi = level_set('fixarea', phi, A0, dx, dy, ep);
        end
    end

    t = t + dt;

    %% ---- 10. diagnostics ---------------------------------------------
    divv = (u(2:end,:)-u(1:end-1,:))/dx + (v(:,2:end)-v(:,1:end-1))/dy;
    if ~isempty(opt.solid), divv(logical(opt.solid)) = 0; end
    A = level_set('area', phi, dx, dy, ep);

    hist.t(n)=t; hist.umax(n)=max(abs([u(:);v(:)]));
    hist.divMax(n)=max(abs(divv(:))); hist.areaErr(n)=(A-A0)/A0;
    band = abs(phi) < 2*dx;
    if any(band(:)), hist.gmag(n) = mean(gmag(band)); else, hist.gmag(n)=NaN; end

    if ~all(isfinite(u(:))) || ~all(isfinite(phi(:)))
        error('ns2phase:diverged','Diverged at step %d (t = %.3e).', n, t);
    end
    if opt.report && (mod(n,opt.report)==0 || n==1 || n==nsteps)
        fprintf('    %6d   %10.4e  %10.3e  %10.3e  %+9.2e   %8.4f\n', ...
                n, t, hist.umax(n), hist.divMax(n), hist.areaErr(n), hist.gmag(n));
    end
end

S.u=u; S.v=v; S.p=p; S.phi=phi; S.G=G; S.dt=dt; S.t=t; S.steps=nsteps;
S.uc = 0.5*(u(1:end-1,:)+u(2:end,:));
S.vc = 0.5*(v(:,1:end-1)+v(:,2:end));
S.hist=hist; S.opt=opt; S.eps=ep;
S.divMax = hist.divMax(end);
S.areaErr = hist.areaErr(end);
S.dtLimit = names{which};

if opt.report
    fprintf('\n    final: area error %+.3e, max|div u| %.2e, binding limit %s\n\n', ...
            S.areaErr, S.divMax, S.dtLimit);
end
end

%% ======================================================= helpers ==========
function t = pbc(b)
if strcmp(b.type,'outflow'), t='D'; else, t='N'; end
end

function [ue, ve] = pad_vel(u, v, opt)
u(1,:)   = nrm(opt.bc.left,  'u', u(2,:));
u(end,:) = nrm(opt.bc.right, 'u', u(end-1,:));
ue = zeros(size(u,1), size(u,2)+2);
ue(:,2:end-1) = u;
ue(:,1)   = tang(opt.bc.bottom,'u',u(:,1));
ue(:,end) = tang(opt.bc.top,   'u',u(:,end));

v(:,1)   = nrm(opt.bc.bottom,'v', v(:,2));
v(:,end) = nrm(opt.bc.top,   'v', v(:,end-1));
ve = zeros(size(v,1)+2, size(v,2));
ve(2:end-1,:) = v;
ve(1,:)   = tang(opt.bc.left, 'v',v(1,:));
ve(end,:) = tang(opt.bc.right,'v',v(end,:));
end

function val = nrm(b, comp, inner)
switch b.type
    case {'wall','moving','slip'}, val = zeros(size(inner));
    case 'inflow'
        % ORIENTATION TRAP: an inflow profile is naturally built as a COLUMN
        % (from G.yu), but the slice it must fill, us(1,:), is a ROW.  Writing
        % profile .* ones(size(inner)) then broadcasts a column against a row
        % and silently produces an ny-by-ny MATRIX, which reshape rejects with
        % "Number of elements must not change".  Reshape the profile itself.
        p = b.(comp);
        if isscalar(p)
            val = p*ones(size(inner));
        elseif numel(p) == numel(inner)
            val = reshape(p, size(inner));
        else
            error('ns2phase:bcSize', ...
                  'Inflow profile has %d elements but the boundary needs %d.', ...
                  numel(p), numel(inner));
        end
    case 'outflow', val = inner;
    otherwise, error('ns2phase:bc','Unknown BC ''%s''.',b.type);
end
end

function g = tang(b, comp, inner)
switch b.type
    case 'wall',   g = -inner;
    case 'moving', if isfield(b,comp), w=b.(comp); else, w=0; end, g = 2*w-inner;
    case {'slip','outflow'}, g = inner;
    case 'inflow', g = -inner;
    otherwise, error('ns2phase:bc','Unknown BC ''%s''.',b.type);
end
end

function [us, vs] = set_normal(us, vs, opt)
us(1,:)   = nrm(opt.bc.left,  'u', us(2,:));
us(end,:) = nrm(opt.bc.right, 'u', us(end-1,:));
vs(:,1)   = nrm(opt.bc.bottom,'v', vs(:,2));
vs(:,end) = nrm(opt.bc.top,   'v', vs(:,end-1));
end
