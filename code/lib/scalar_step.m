function c = scalar_step(c, u, v, D, dt, dx, dy, lim, solid)
%SCALAR_STEP  One explicit step of 2-D advection-diffusion on a MAC grid.
%
%   c = scalar_step(c, u, v, D, dt, dx, dy, lim, solid)
%
%   c      nx-by-ny cell-centred concentration
%   u, v   staggered face velocities, (nx+1)xny and nx x(ny+1)
%   lim    flux limiter name (see LIMITER); 'vanleer' is the right default
%   solid  optional nx-by-ny logical mask; no flux crosses a solid face
%
%   CONSERVATIVE FLUX FORM.  Advective fluxes are evaluated at the SAME faces
%   where the velocities live -- that is the payoff of the staggered grid, and
%   it means total scalar mass is conserved to machine precision (no
%   interpolation of velocity is required anywhere).
%
%   Walls and solid faces are zero-flux for both advection and diffusion, so
%   sum(c) changes only through open inflow/outflow boundaries.
%
%   WHY THE LIMITER MATTERS HERE MORE THAN ANYWHERE ELSE
%   At the Pe ~ 200 of a microchannel, first-order upwind produces numerical
%   diffusion COMPARABLE TO THE PHYSICAL DIFFUSIVITY (Tutorial 5).  A mixer
%   simulated that way appears to work beautifully and is entirely wrong.
%   vanleer costs the same per step and produces ~750x less.
%
%   See also LIMITER, MIXING_STUDY, TAYLOR_DISPERSION, NS_SOLVER.

if nargin < 8 || isempty(lim), lim = 'vanleer'; end
if nargin < 9, solid = []; end

[nx, ny] = size(c);
hasSolid = ~isempty(solid);
if hasSolid, solid = logical(solid); end

%% -------------------------------------------------- advective x-fluxes ---
% Face i (i = 2..nx) sits between cells i-1 and i.  Build the TVD flux there.
uF = u(2:nx, :);                                  % (nx-1) x ny interior faces
cL = c(1:nx-1, :);   cR = c(2:nx, :);             % cells either side
dcF = cR - cL;

% Gradient ratio r, using the cell UPWIND of the face.
cLL = [c(1,:); c(1:nx-2,:)];                      % c_{i-2}, clamped at the edge
cRR = [c(3:nx,:); c(nx,:)];                       % c_{i+1}, clamped
rPos = (cL - cLL) ./ dcF;
rNeg = (cRR - cR) ./ dcF;

Cno = abs(uF)*dt/dx;
pos = uF >= 0;
r   = rPos.*pos + rNeg.*(~pos);
psi = limiter(r, lim);

cUp = cL.*pos + cR.*(~pos);
Fx  = uF.*cUp + 0.5*abs(uF).*(1 - Cno).*psi.*dcF;

% Diffusive x-flux, same faces
Fx = Fx - D*dcF/dx;

%% -------------------------------------------------- advective y-fluxes ---
vF  = v(:, 2:ny);
cB  = c(:, 1:ny-1);  cT = c(:, 2:ny);
dcG = cT - cB;

cBB = [c(:,1), c(:,1:ny-2)];
cTT = [c(:,3:ny), c(:,ny)];
rPos = (cB - cBB) ./ dcG;
rNeg = (cTT - cT) ./ dcG;

Cnv = abs(vF)*dt/dy;
pos = vF >= 0;
r   = rPos.*pos + rNeg.*(~pos);
psi = limiter(r, lim);

cUp = cB.*pos + cT.*(~pos);
Fy  = vF.*cUp + 0.5*abs(vF).*(1 - Cnv).*psi.*dcG;
Fy  = Fy - D*dcG/dy;

%% ------------------------------------------------------- block solid -----
if hasSolid
    Fx(solid(1:nx-1,:) | solid(2:nx,:)) = 0;
    Fy(solid(:,1:ny-1) | solid(:,2:ny)) = 0;
end

%% ------------------------------------------------------------ update -----
% Pad with zero flux at the outer walls; open boundaries are handled by the
% caller (which overwrites inlet cells and extrapolates at the outlet).
FxA = [zeros(1,ny); Fx; zeros(1,ny)];
FyA = [zeros(nx,1), Fy, zeros(nx,1)];

c = c - (dt/dx)*(FxA(2:end,:) - FxA(1:end-1,:)) ...
      - (dt/dy)*(FyA(:,2:end) - FyA(:,1:end-1));

if hasSolid, c(solid) = 0; end
end
