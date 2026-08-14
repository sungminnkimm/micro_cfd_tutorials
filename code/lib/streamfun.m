function psi = streamfun(S)
%STREAMFUN  Streamfunction from a staggered velocity field.
%
%   psi = streamfun(S)   S from NS_SOLVER
%
%   For 2-D incompressible flow, u = dpsi/dy and v = -dpsi/dx, so psi can be
%   recovered by integrating u along y.  Returned at cell centres (nx x ny).
%
%   Streamlines are contours of psi, and because the flow is divergence-free
%   those contours are exact particle paths in steady flow -- much more
%   informative than a quiver plot for spotting recirculation zones.
%
%   The additive constant is fixed by psi = 0 at the bottom wall.  Integrating
%   u along y is well-posed ONLY because div(u) = 0 to machine precision on the
%   staggered grid; on a colocated grid the result would depend on the path you
%   integrated along, which is a good way to notice you have a divergence
%   problem.
%
%   See also NS_SOLVER, CAVITY.

uc = S.uc;                 % nx x ny, cell-centred u
dy = S.G.dy;

psi = zeros(size(uc));
psi(:,1) = uc(:,1)*dy/2;                    % from the wall to the first centre
for j = 2:size(uc,2)
    psi(:,j) = psi(:,j-1) + 0.5*(uc(:,j-1) + uc(:,j))*dy;   % trapezoid
end
end
