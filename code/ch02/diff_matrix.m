function [D1, D2] = diff_matrix(N, h, bc)
%DIFF_MATRIX  Sparse 1-D first and second derivative operators (cell-centred).
%
%   [D1,D2] = diff_matrix(N,h,bc) returns N-by-N sparse matrices acting on a
%   column vector of N CELL-CENTRED values with spacing h.  Cell j is centred
%   at x = (j-1/2)*h, so the domain boundaries sit half a cell outside the
%   first and last unknowns.
%
%   bc = 'dirichlet0'  u = 0 at both walls      -> ghost u_0 = -u_1
%      = 'neumann0'    du/dn = 0 at both walls  -> ghost u_0 = +u_1
%      = 'periodic'    wraps around
%
%   The ghost-node relations are enforced by folding the ghost coefficient back
%   into the first/last diagonal entry.  Derivations (left wall, ghost u_0):
%
%     dirichlet0, u_0 = -u_1:
%        D2 row 1 = (u_0 - 2u_1 + u_2)/h^2 = (-3u_1 + u_2)/h^2
%        D1 row 1 = (u_2 - u_0)/(2h)       = (+u_1 + u_2)/(2h)
%     neumann0, u_0 = +u_1:
%        D2 row 1 = (u_0 - 2u_1 + u_2)/h^2 = ( -u_1 + u_2)/h^2
%        D1 row 1 = (u_2 - u_0)/(2h)       = (-u_1 + u_2)/(2h)
%
%   Accuracy note: D2 with these ghost relations is globally second order --
%   the standard finite-volume result, and the reason cell-centring is worth
%   the half-cell bookkeeping.  The boundary ROWS of D1 are only first order;
%   this course never differentiates a scalar right at a wall, so it does not
%   bite, but do not reuse D1's corner rows in a high-accuracy setting.
%
%   Example:
%       [D1,D2] = diff_matrix(8,0.1,'dirichlet0');
%       full(D2), spy(D2)
%
%   See also CONVERGENCE_STUDY, POISSON2D.

e  = ones(N,1);
D1 = spdiags([-e  0*e  e], [-1 0 1], N, N) / (2*h);
D2 = spdiags([ e -2*e  e], [-1 0 1], N, N) / h^2;

switch lower(bc)
    case 'periodic'
        D1(1,N) = -1/(2*h);   D1(N,1) =  1/(2*h);
        D2(1,N) =  1/h^2;     D2(N,1) =  1/h^2;

    case 'dirichlet0'
        D2(1,1) = -3/h^2;     D2(N,N) = -3/h^2;
        D1(1,1) =  1/(2*h);   D1(N,N) = -1/(2*h);

    case 'neumann0'
        D2(1,1) = -1/h^2;     D2(N,N) = -1/h^2;
        D1(1,1) = -1/(2*h);   D1(N,N) =  1/(2*h);

    otherwise
        error('diff_matrix:bc','Unknown bc ''%s''. Use dirichlet0|neumann0|periodic.', bc);
end
end
