function psi = limiter(r, name)
%LIMITER  Flux limiter functions psi(r) for TVD advection schemes.
%
%   psi = limiter(r, name)
%
%   r is the ratio of consecutive gradients, r_i = (c_i - c_{i-1})/(c_{i+1} - c_i),
%   a local smoothness sensor: r ~ 1 means smooth, r <= 0 means a local extremum.
%
%   name = 'upwind'    psi = 0                    most diffusive, 1st order
%        = 'lw'        psi = 1                    Lax-Wendroff, 2nd order, wiggly
%        = 'vanleer'   (r+|r|)/(1+|r|)            smooth, good default
%        = 'superbee'  max(0,min(2r,1),min(r,2))  sharpest, can over-steepen
%        = 'minmod'    max(0,min(r,1))            most cautious
%
%   Every TVD limiter satisfies psi(r) = 0 for r <= 0 (kill the correction at
%   extrema, so no new maxima can be created -- this is what prevents wiggles)
%   and psi(1) = 1 (full second-order accuracy where the solution is smooth).
%
%   NaN GUARD: r is Inf or NaN wherever c_{i+1} == c_i (a flat region divides
%   by zero).  Those must map to psi = 0.  An unguarded NaN contaminates the
%   entire field within one timestep and is the classic first-run failure of a
%   hand-written TVD scheme -- if your solution goes all-NaN in one step rather
%   than growing over many, look here before you look at stability.
%
%   See also ADVECT1D, ADVECT_TVD.

switch lower(name)
    case 'upwind',   psi = zeros(size(r));
    case 'lw',       psi = ones(size(r));
    case 'vanleer',  psi = (r + abs(r))./(1 + abs(r));
    case 'superbee', psi = max(0, max(min(2*r,1), min(r,2)));
    case 'minmod',   psi = max(0, min(r,1));
    otherwise
        error('limiter:name', ...
              'Unknown limiter ''%s''. Use upwind|lw|vanleer|superbee|minmod.', name);
end

psi(~isfinite(psi)) = 0;
end
