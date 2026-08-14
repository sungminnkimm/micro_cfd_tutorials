function cNew = advect_tvd(c, u, dt, dx, lim)
%ADVECT_TVD  One explicit TVD step of 1-D advection on a periodic domain.
%
%   cNew = advect_tvd(c, u, dt, dx, lim)
%
%   c    column vector of cell-averaged values
%   u    scalar advection velocity (either sign)
%   lim  limiter name, see LIMITER
%
%   CONSERVATIVE FLUX FORM
%   ----------------------
%       c_i^{n+1} = c_i^n - (dt/dx) * (F_{i+1/2} - F_{i-1/2})
%
%   Because every interior face flux appears once with each sign, the total
%   sum(c) is conserved to MACHINE PRECISION regardless of the limiter.  That
%   is worth having for its own sake: it turns "is my scalar leaking?" into a
%   question with a crisp yes/no answer instead of a judgement call.
%
%   The face flux is upwind plus a limited anti-diffusive correction:
%
%       F_{i+1/2} = u*c_upwind + 0.5*|u|*(1 - |C|)*psi(r_i)*(c_{i+1} - c_i)
%
%   The (1 - |C|) factor is not cosmetic.  It vanishes at |C| = 1, where the
%   scheme becomes EXACT (a parcel moves precisely one cell per step and the
%   update is a pure shift).  It is also why a larger Courant number gives LESS
%   numerical diffusion for upwind -- see Exercise 5.2, which is worth doing
%   because the result is genuinely counter-intuitive.
%
%   See also LIMITER, ADVECT1D, NUM_DIFFUSION.

c   = c(:);
Cno = u*dt/dx;

cm1 = circshift(c, 1);      % c_{i-1}
cp1 = circshift(c,-1);      % c_{i+1}
cp2 = circshift(c,-2);      % c_{i+2}

dcp = cp1 - c;              % gradient across face i+1/2

if u >= 0
    r      = (c - cm1) ./ dcp;      % upwind side is to the left
    cUpw   = c;
else
    r      = (cp2 - cp1) ./ dcp;    % upwind side is to the right
    cUpw   = cp1;
end

psi = limiter(r, lim);

F = u.*cUpw + 0.5*abs(u)*(1 - abs(Cno)).*psi.*dcp;

cNew = c - (dt/dx)*(F - circshift(F,1));
end
