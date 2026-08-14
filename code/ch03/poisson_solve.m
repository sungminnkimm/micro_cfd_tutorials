function phi = poisson_solve(S, f, bcVals)
%POISSON_SOLVE  Solve del^2(phi) = f using a factorization from POISSON2D.
%
%   phi = poisson_solve(S, f, bcVals)
%
%   f       nx-by-ny array, the source term.
%   bcVals  1x4 cell {xlo,xhi,ylo,yhi}.  For a 'D' wall the entry is the wall
%           VALUE phi_w; for an 'N' wall it is the OUTWARD normal derivative
%           g = dphi/dn.  Each entry may be a scalar or a vector along the wall
%           (length ny for x-walls, nx for y-walls).
%
%   Boundary contributions move to the right-hand side.  Derivations, left wall
%   (ghost cell 0 sits half a cell outside cell 1):
%
%     Dirichlet, phi_0 = 2*phi_w - phi_1:
%        (phi_0 - 2phi_1 + phi_2)/dx^2 = (-3phi_1 + phi_2)/dx^2 + 2*phi_w/dx^2
%        => b = f - 2*phi_w/dx^2
%
%     Neumann with OUTWARD normal (-x at the left wall), dphi/dn = g:
%        dphi/dx|wall = -g  =>  phi_0 = phi_1 + g*dx
%        (phi_0 - 2phi_1 + phi_2)/dx^2 = (-phi_1 + phi_2)/dx^2 + g/dx
%        => b = f - g/dx
%
%   Both walls give a MINUS sign because g is referred to the outward normal.
%   That symmetry is not cosmetic: summing the equations over all cells then
%   gives exactly the compatibility condition  int(f)dV = oint(g)dS.  If your
%   signs are inconsistent, an all-Neumann problem will silently violate it.
%
%   See also POISSON2D.

b = f;

if S.bcTypes{1}=='D', b(1,:)   = b(1,:)   - 2*bcVals{1}(:).'/S.dx^2;
else,                 b(1,:)   = b(1,:)   -   bcVals{1}(:).'/S.dx;   end

if S.bcTypes{2}=='D', b(end,:) = b(end,:) - 2*bcVals{2}(:).'/S.dx^2;
else,                 b(end,:) = b(end,:) -   bcVals{2}(:).'/S.dx;   end

if S.bcTypes{3}=='D', b(:,1)   = b(:,1)   - 2*bcVals{3}(:)/S.dy^2;
else,                 b(:,1)   = b(:,1)   -   bcVals{3}(:)/S.dy;     end

if S.bcTypes{4}=='D', b(:,end) = b(:,end) - 2*bcVals{4}(:)/S.dy^2;
else,                 b(:,end) = b(:,end) -   bcVals{4}(:)/S.dy;     end

b = b(:);
if isfield(S,'mask') && any(S.mask(:))
    b(S.mask(:)) = 0;               % solid rows are the identity
end
if S.pinned, b(S.pinIdx) = 0; end   % must match the pinned row in A

phi = S.Q * (S.U \ (S.L \ (S.P * (S.R \ b))));
phi = reshape(phi, S.nx, S.ny);
end
