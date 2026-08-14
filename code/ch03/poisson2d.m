function S = poisson2d(nx, ny, dx, dy, bcTypes, mask)
%POISSON2D  Assemble and factorize the 2-D cell-centred Laplacian.
%
%   S = poisson2d(nx,ny,dx,dy,bcTypes)          plain rectangular domain
%   S = poisson2d(nx,ny,dx,dy,bcTypes,mask)     with solid cells excluded
%
%   Solve with:   phi = poisson_solve(S, f, bcVals);
%
%   bcTypes is a 1x4 cell {xlo, xhi, ylo, yhi}, each 'D' (Dirichlet) or 'N'
%   (Neumann).  mask is an optional nx-by-ny logical, TRUE in solid cells.
%
%   Unknowns are flattened column-major, k = i + (j-1)*nx, which is what
%   phi(:) / reshape(phi,nx,ny) do natively.  Index order is (x,y), so build
%   coordinate arrays with NDGRID, not meshgrid, and transpose before plotting.
%
%   ASSEMBLY BY FLUX BALANCE
%   ------------------------
%   Each row is the sum over the cell's four faces of (phi_nb - phi_c)*c, with
%   c = 1/dx^2 or 1/dy^2.  Every boundary type is then just a rule for one face:
%
%     interior face      A(k,knb) += c ;  A(k,k) -= c
%     Dirichlet wall     A(k,k)   -= 2c        (ghost = 2*phi_w - phi_c)
%     Neumann wall       nothing               (zero flux through the face)
%     SOLID neighbour    nothing               (also zero flux)
%
%   Note the last two lines are identical: a solid face IS a Neumann face.
%   That is the whole trick, and it is why masking has to happen HERE rather
%   than by zeroing velocities afterwards.  If the Poisson operator does not
%   know about the body, the projection cannot make the flow divergence-free
%   around it -- you get div ~ 1e2 instead of 1e-14 (Tutorial 8.3).
%
%   Assembly is expensive; the solve is cheap.  Call poisson2d ONCE outside
%   your time loop and poisson_solve every step.
%
%   See also POISSON_SOLVE, POISSON_MMS, NS_SOLVER, MAKE_MASK.

if numel(bcTypes)~=4
    error('poisson2d:bc','bcTypes must be a 1x4 cell {xlo,xhi,ylo,yhi}.');
end
if nargin < 6 || isempty(mask)
    mask = false(nx,ny);
elseif ~isequal(size(mask),[nx ny])
    error('poisson2d:mask','mask must be %dx%d.', nx, ny);
end
mask = logical(mask);

cx = 1/dx^2;  cy = 1/dy^2;
N  = nx*ny;

% Triplet lists.  5 entries per cell is the upper bound.
I = zeros(5*N,1);  J = I;  V = I;  ptr = 0;
    function add(r,c,v)
        ptr = ptr+1;  I(ptr)=r;  J(ptr)=c;  V(ptr)=v;
    end

for j = 1:ny
    for i = 1:nx
        k = i + (j-1)*nx;

        if mask(i,j)
            add(k,k,1);              % solid cell: phi = 0, decoupled
            continue
        end

        d = 0;

        % ---- x-low face
        if i == 1
            if bcTypes{1}=='D', d = d - 2*cx; end          % 'N' adds nothing
        elseif mask(i-1,j)
            % solid neighbour: no flux, nothing to add
        else
            add(k, k-1, cx);  d = d - cx;
        end

        % ---- x-high face
        if i == nx
            if bcTypes{2}=='D', d = d - 2*cx; end
        elseif mask(i+1,j)
        else
            add(k, k+1, cx);  d = d - cx;
        end

        % ---- y-low face
        if j == 1
            if bcTypes{3}=='D', d = d - 2*cy; end
        elseif mask(i,j-1)
        else
            add(k, k-nx, cy);  d = d - cy;
        end

        % ---- y-high face
        if j == ny
            if bcTypes{4}=='D', d = d - 2*cy; end
        elseif mask(i,j+1)
        else
            add(k, k+nx, cy);  d = d - cy;
        end

        if d == 0
            % Fully enclosed fluid cell with no Dirichlet anywhere touching it
            % (e.g. a 1-cell pocket).  Leave it decoupled rather than singular.
            d = 1;
        end
        add(k,k,d);
    end
end

A = sparse(I(1:ptr), J(1:ptr), V(1:ptr), N, N);

% All-Neumann and no solid Dirichlet anywhere -> the constant is a null mode.
% Pin the first FLUID cell.  Legitimate: only grad(phi) is ever used.
S.pinned = all(strcmp(bcTypes,'N'));
if S.pinned
    kp = find(~mask(:), 1);
    A(kp,:) = 0;  A(kp,kp) = 1;
    S.pinIdx = kp;
else
    S.pinIdx = [];
end

S.nx = nx;  S.ny = ny;  S.dx = dx;  S.dy = dy;
S.bcTypes = bcTypes;  S.mask = mask;
S.A = A;
[S.L, S.U, S.P, S.Q, S.R] = lu(A);      % factorize ONCE
end
