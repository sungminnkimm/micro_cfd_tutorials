function cp = padarray_neumann(c)
%PADARRAY_NEUMANN  Add a one-cell ghost ring enforcing zero normal flux.
%
%   cp = padarray_neumann(c) returns an (nx+2)-by-(ny+2) array whose interior
%   cp(2:end-1,2:end-1) is c and whose ghost ring mirrors the adjacent interior
%   value, i.e. c_ghost = c_adjacent.  This makes the one-sided difference
%   across the wall vanish, which is exactly du/dn = 0.
%
%   Written as an explicit function because ghost-cell handling is where most
%   finite-difference bugs live, and having ONE place that does it means there
%   is one place to check.  Corner ghosts are never used by the five-point
%   stencil, but are filled consistently so the array is safe to differentiate
%   twice.
%
%   Equivalent to padarray(c,[1 1],'replicate') from the Image Processing
%   Toolbox; written out here so the course needs no toolboxes.
%
%   See also DIFFUSION2D.

cp = zeros(size(c,1)+2, size(c,2)+2, 'like', c);
cp(2:end-1, 2:end-1) = c;

cp(1,   2:end-1) = c(1,  :);      % x-low  wall
cp(end, 2:end-1) = c(end,:);      % x-high wall
cp(2:end-1, 1  ) = c(:,  1);      % y-low  wall
cp(2:end-1, end) = c(:,end);      % y-high wall

cp(1,1)     = c(1,1);             % corners (unused by 5-point stencils)
cp(1,end)   = c(1,end);
cp(end,1)   = c(end,1);
cp(end,end) = c(end,end);
end
