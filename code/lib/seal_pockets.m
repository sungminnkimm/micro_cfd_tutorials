function [maskOut, nSealed] = seal_pockets(mask, seed)
%SEAL_POCKETS  Fill fluid regions that are not connected to the main flow.
%
%   [maskOut, nSealed] = seal_pockets(mask, seed)
%
%   mask    nx-by-ny logical, TRUE = solid
%   seed    optional [i j] fluid cell to flood from; default is the first fluid
%           cell in column 1 (i.e. the inlet)
%
%   Returns the mask with every fluid cell NOT reachable from the seed turned
%   solid, plus a count of how many were sealed.
%
%   WHY THIS MATTERS
%   ----------------
%   A staircased diagonal channel (see MAKE_MASK) can leave small pockets of
%   fluid completely surrounded by solid.  Such a pocket is an ISOLATED
%   all-Neumann sub-domain in the pressure Poisson equation: the constant is a
%   null mode there, the pinning cell is somewhere else entirely, and the solve
%   returns garbage pressure inside the pocket.  The garbage feeds back through
%   the corrector and the whole simulation diverges -- typically after a few
%   hundred steps, which makes it look like a timestep problem when it is not.
%
%   The symptom is distinctive: "diverged" at a dt COMFORTABLY INSIDE the
%   stability limits.  If your dt is half the viscous limit and it still blows
%   up, stop tuning the timestep and look at your geometry.
%
%   Uses a breadth-first flood fill, so no Image Processing Toolbox is needed.
%
%   Example:
%       solid = make_mask(G, {'below',@(x) f(x)-w}, {'above',@(x) f(x)+w});
%       [solid, n] = seal_pockets(solid);
%       if n, fprintf('sealed %d isolated cells\n', n); end
%
%   See also MAKE_MASK, NS_SOLVER, POISSON2D.

mask = logical(mask);
[nx, ny] = size(mask);

if nargin < 2 || isempty(seed)
    j = find(~mask(1,:), 1);
    if isempty(j)
        % Inlet column fully blocked; fall back to any fluid cell.
        k = find(~mask(:), 1);
        if isempty(k)
            maskOut = mask;  nSealed = 0;  return
        end
        [i, j] = ind2sub([nx ny], k);
        seed = [i j];
    else
        seed = [1 j];
    end
end

reach = false(nx, ny);
if mask(seed(1), seed(2))
    error('seal_pockets:seed','Seed cell (%d,%d) is solid.', seed(1), seed(2));
end

% Breadth-first flood fill over 4-connected fluid cells.
queue = zeros(nx*ny, 2);
queue(1,:) = seed;  head = 1;  tail = 1;
reach(seed(1), seed(2)) = true;

while head <= tail
    i = queue(head,1);  j = queue(head,2);  head = head + 1;
    nb = [i-1 j; i+1 j; i j-1; i j+1];
    for q = 1:4
        a = nb(q,1);  b = nb(q,2);
        if a >= 1 && a <= nx && b >= 1 && b <= ny && ~mask(a,b) && ~reach(a,b)
            reach(a,b) = true;
            tail = tail + 1;
            queue(tail,:) = [a b];
        end
    end
end

isolated = ~mask & ~reach;
nSealed  = sum(isolated(:));
maskOut  = mask | isolated;
end
