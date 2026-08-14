function [k, nx_, ny_, gmag] = curvature(phi, dx, dy)
%CURVATURE  Interface curvature and normal from a level-set function.
%
%   [k, nx, ny, gmag] = curvature(phi, dx, dy)
%
%   k     total curvature, div(grad(phi)/|grad(phi)|)
%   nx,ny unit normal components
%   gmag  |grad(phi)| -- should be ~1 for a proper signed distance function,
%         and watching it drift away from 1 is how you know reinitialization is
%         overdue (Tutorial 11.3)
%
%   In 2-D the expanded form avoids differentiating the normalized gradient
%   twice, which is noisier:
%
%       k = (phi_xx*phi_y^2 - 2*phi_x*phi_y*phi_xy + phi_yy*phi_x^2)
%           / (phi_x^2 + phi_y^2)^(3/2)
%
%   SIGN CONVENTION -- read this, it is easy to get backwards.
%
%   With phi > 0 INSIDE a droplet, grad(phi) points INWARD, so the raw formula
%   above returns -1/R for a circle:  for phi = R - r we get grad(phi) = -rhat,
%   and div(-rhat) = -1/r in 2-D.  That is geometrically correct but the
%   opposite of what "the curvature of a droplet" usually means, and using it
%   directly in the CSF force makes surface tension push the droplet APART.
%
%   So this function NEGATES it and returns k = +1/R for a droplet, which:
%     * matches Young-Laplace as usually written, dp = sigma*k, with the higher
%       pressure inside;
%     * makes the CSF force  sigma*k*delta*grad(phi)  point INWARD, since
%       grad(phi) already points inward.  The two sign conventions cancel.
%
%   In 2-D this means dp = sigma/R, NOT 2*sigma/R -- a 2-D droplet is an
%   infinite cylinder, not a sphere.  Comparing a 2-D result against the
%   spherical formula is one of the most common factor-of-2 errors in two-phase
%   CFD (Tutorial 10.2).  Getting the SIGN wrong is a more spectacular one: the
%   symptom is a negative pressure jump and a droplet that tears itself apart.
%
%   ACCURACY WARNING: curvature is a SECOND derivative, so it amplifies noise in
%   phi by 1/dx^2.  This is why phi must be kept smooth (a signed distance
%   function), and it is the root cause of spurious currents.  Use >= 10 cells
%   per droplet radius, 20+ if you care about the pressure jump.
%
%   Boundaries use replicated (zero-gradient) ghost cells.
%
%   See also LEVEL_SET, YOUNG_LAPLACE, NS2PHASE.

p = padarray_neumann(phi);

% First derivatives (central, on the padded array -> same size as phi)
px = (p(3:end,2:end-1) - p(1:end-2,2:end-1)) / (2*dx);
py = (p(2:end-1,3:end) - p(2:end-1,1:end-2)) / (2*dy);

% Second derivatives
pxx = (p(3:end,2:end-1) - 2*phi + p(1:end-2,2:end-1)) / dx^2;
pyy = (p(2:end-1,3:end) - 2*phi + p(2:end-1,1:end-2)) / dy^2;

% Mixed derivative needs the corners of the padded array
pxy = (p(3:end,3:end) - p(3:end,1:end-2) ...
     - p(1:end-2,3:end) + p(1:end-2,1:end-2)) / (4*dx*dy);

gmag = sqrt(px.^2 + py.^2);
den  = max(gmag.^3, eps);

% Negated so that phi > 0 inside a droplet gives k = +1/R -- see the sign
% discussion in the header.
k = -(pxx.*py.^2 - 2*px.*py.*pxy + pyy.*px.^2) ./ den;

g    = max(gmag, eps);
nx_  = px ./ g;
ny_  = py ./ g;
end
