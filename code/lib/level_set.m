function out = level_set(action, varargin)
%LEVEL_SET  Level-set interface tracking utilities.
%
%   phi = level_set('circle', G, xc, yc, R)     signed distance to a circle
%   phi = level_set('advect', phi, u, v, dt, dx, dy)     one TVD advection step
%   phi = level_set('reinit', phi, dx, dy, nIter)        restore |grad phi| = 1
%   H   = level_set('heaviside', phi, eps)      smoothed Heaviside
%   d   = level_set('delta',     phi, eps)      smoothed delta
%   A   = level_set('area',      phi, dx, dy, eps)   enclosed area (phi > 0)
%   level_set('selftest')                       Zalesak disk rotation test
%
%   CONVENTION: phi > 0 INSIDE the droplet, phi < 0 outside, phi = 0 on the
%   interface.  Signed distance means |grad phi| = 1, which is what makes
%   normals and curvature clean -- and keeping that property is the whole
%   difficulty (Tutorial 11.3).
%
%   REINITIALIZATION solves, in pseudo-time,
%       dphi/dtau = sign(phi0) * (1 - |grad phi|)
%   which drives |grad phi| -> 1 while holding the zero contour fixed.  It uses
%   a SMOOTHED sign function, because the discrete sign is not exactly zero at
%   the interface and the resulting drift is the main source of level-set mass
%   loss.  Reinitialize too often and you lose mass; too rarely and curvature
%   is corrupted.  There is no universally right setting -- tune it against the
%   area diagnostic.
%
%   See also CURVATURE, NS2PHASE, STATIC_DROPLET.

switch lower(action)
    case 'circle'
        [G, xc, yc, R] = varargin{1:4};
        out = R - sqrt((G.Xp-xc).^2 + (G.Yp-yc).^2);

    case 'advect'
        out = ls_advect(varargin{:});

    case 'reinit'
        out = ls_reinit(varargin{:});

    case 'heaviside'
        [phi, ep] = varargin{1:2};
        out = 0.5*(1 + phi/ep + sin(pi*phi/ep)/pi);
        out(phi < -ep) = 0;
        out(phi >  ep) = 1;

    case 'delta'
        [phi, ep] = varargin{1:2};
        out = (1 + cos(pi*phi/ep))/(2*ep);
        out(abs(phi) > ep) = 0;

    case 'area'
        [phi, dx, dy, ep] = varargin{1:4};
        H = level_set('heaviside', phi, ep);
        out = sum(H(:))*dx*dy;

    case 'fixarea'
        % Global area correction: shift phi by a constant so the enclosed area
        % returns to its target.  Reinitialization drifts the zero contour
        % (the discrete sign function is not exactly zero at the interface),
        % and at ~0.1% per call over several hundred calls that is enough to
        % dissolve a droplet completely.  This is the standard cheap remedy.
        %
        % Newton step: dA/d(shift) = integral of delta(phi) = the perimeter.
        %
        % HONEST LIMITATION: this fixes TOTAL volume, not its local
        % distribution, and it shifts the interface uniformly rather than where
        % the loss actually occurred.  It is a patch on a known weakness of
        % level sets, not a cure.  If exact droplet volumes are your
        % quantity of interest, use VOF or CLSVOF instead (Tutorial 11.1).
        [phi, Atarget, dx, dy, ep] = varargin{1:5};
        for it = 1:4
            H = level_set('heaviside', phi, ep);
            A = sum(H(:))*dx*dy;
            d = level_set('delta', phi, ep);
            Per = sum(d(:))*dx*dy;
            if Per < eps || abs(A-Atarget) < 1e-12*max(Atarget,eps), break; end
            phi = phi + (Atarget - A)/Per;
        end
        out = phi;

    case 'selftest'
        out = zalesak_test(varargin{:});

    otherwise
        error('level_set:action','Unknown action ''%s''.', action);
end
end

%% ==========================================================================
function phi = ls_advect(phi, u, v, dt, dx, dy)
%LS_ADVECT  One TVD (van Leer) advection step of phi on a MAC velocity field.
%   phi is smooth, so the wiggle problem of Tutorial 5 is absent -- but
%   numerical diffusion here directly MOVES THE INTERFACE, which is worse.
%   Reuse the conservative scalar step; for a divergence-free velocity field
%   the conservative and advective forms coincide.
phi = scalar_step(phi, u, v, 0, dt, dx, dy, 'vanleer');
end

%% ==========================================================================
function phi = ls_reinit(phi, dx, dy, nIter)
%LS_REINIT  Restore the signed-distance property.
if nargin < 4, nIter = 5; end

phi0 = phi;
% Smoothed sign: phi0 / sqrt(phi0^2 + |grad phi0|^2 dx^2).  Reduces (but does
% not eliminate) the zero-contour drift that costs mass.
p0 = padarray_neumann(phi0);
gx = (p0(3:end,2:end-1)-p0(1:end-2,2:end-1))/(2*dx);
gy = (p0(2:end-1,3:end)-p0(2:end-1,1:end-2))/(2*dy);
S  = phi0 ./ sqrt(phi0.^2 + (gx.^2+gy.^2)*dx^2 + eps);

dtau = 0.5*min(dx,dy);

for it = 1:nIter
    p = padarray_neumann(phi);

    % One-sided differences for a Godunov upwind |grad phi|
    a = (phi - p(1:end-2,2:end-1))/dx;      % backward x
    b = (p(3:end,2:end-1) - phi)/dx;        % forward  x
    c = (phi - p(2:end-1,1:end-2))/dy;      % backward y
    d = (p(2:end-1,3:end) - phi)/dy;        % forward  y

    ap = max(a,0); am = min(a,0);
    bp = max(b,0); bm = min(b,0);
    cp = max(c,0); cm = min(c,0);
    dp_ = max(d,0); dm = min(d,0);

    gP = sqrt(max(ap.^2, bm.^2) + max(cp.^2, dm.^2));   % where S > 0
    gM = sqrt(max(am.^2, bp.^2) + max(cm.^2, dp_.^2));  % where S < 0

    G = (S>0).*gP + (S<0).*gM;
    phi = phi - dtau*S.*(G - 1);
end
end

%% ==========================================================================
function res = zalesak_test(N, doReinit)
%ZALESAK_TEST  Slotted disk in solid-body rotation: the standard interface test.
%   After one full revolution the shape should return to its initial state.
%   Isolates interface-tracking error -- no surface tension, no coupling.
if nargin < 1 || isempty(N), N = 100; end
if nargin < 2, doReinit = true; end

L = 1;  G = mac_grid(N,N,L,L);
dx = G.dx;  dy = G.dy;

% Slotted disk centred at (0.5,0.75), radius 0.15, slot 0.05 wide.
% Build a PROPER signed distance function: phi = min(dist_to_disk, -dist_to_slot),
% positive inside the disk and outside the slot.  Hacking the slot in by
% overwriting values (as a first attempt naturally does) leaves phi far from a
% distance function there, and reinitialization then immediately deforms the
% slot -- so the test would be measuring the bad initial condition rather than
% the advection scheme.
R = 0.15;  cx = 0.5;  cy = 0.75;
dDisk = R - sqrt((G.Xp-cx).^2 + (G.Yp-cy).^2);

sw = 0.025;  sTop = cy + 0.09;                    % slot half-width and top
dRect = rect_sdf(G.Xp, G.Yp, cx-sw, cx+sw, cy-2*R, sTop);   % >0 inside slot

phi  = min(dDisk, -dRect);
phi0 = phi;

% Solid-body rotation about the domain centre, omega = 1
om = 1;
u = -om*(G.Yu - L/2);
v =  om*(G.Xv - L/2);

T  = 2*pi/om;
dt = 0.4*dx/max([max(abs(u(:))), max(abs(v(:)))]);
nst = ceil(T/dt);  dt = T/nst;

ep = 1.5*dx;
A0 = level_set('area', phi, dx, dy, ep);

for n = 1:nst
    phi = ls_advect(phi, u, v, dt, dx, dy);
    if doReinit && mod(n,10)==0
        phi = ls_reinit(phi, dx, dy, 3);
    end
end

A1 = level_set('area', phi, dx, dy, ep);
H0 = level_set('heaviside', phi0, ep);
H1 = level_set('heaviside', phi,  ep);

res.shapeErr = sum(abs(H1(:)-H0(:)))*dx*dy / (sum(H0(:))*dx*dy);
res.areaErr  = (A1-A0)/A0;
res.N = N;  res.steps = nst;  res.reinit = doReinit;
res.phi = phi;  res.phi0 = phi0;  res.G = G;
end

%% ==========================================================================
function d = rect_sdf(X, Y, x0, x1, y0, y1)
%RECT_SDF  Signed distance to an axis-aligned rectangle, POSITIVE inside.
dxi = min(X - x0, x1 - X);          % >0 inside the x-slab
dyi = min(Y - y0, y1 - Y);          % >0 inside the y-slab
inside = dxi > 0 & dyi > 0;
d = zeros(size(X));
d(inside) = min(dxi(inside), dyi(inside));
o = ~inside;
d(o) = -sqrt(max(-dxi(o),0).^2 + max(-dyi(o),0).^2);
end
