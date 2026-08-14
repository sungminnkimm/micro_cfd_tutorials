%YOUNG_LAPLACE  Verify the curvature routine, and the pressure it implies.
%
%   Tutorial 10.  Run:  young_laplace
%
%   Step 2 of the Tutorial 11.6 validation ladder.  A circle has k = 1/R
%   EXACTLY, so this isolates curvature error with nothing else in the way.
%   Do it before trusting any two-phase result.
%
%   Also demonstrates the noise sensitivity that motivates reinitialization:
%   curvature is a second derivative, so it amplifies errors in phi by 1/dx^2.
%
%   See also CURVATURE, LEVEL_SET, STATIC_DROPLET.

if ~exist('QUIET','var'), QUIET = false; end

R  = 25e-6;                       % droplet radius
L  = 100e-6;                      % domain
sigma = 0.072;

%% ------------------------------------------- convergence of curvature ----
Ns  = [32 64 128 256];
err = zeros(size(Ns));
fprintf('\n  Curvature of a circle, R = %g um.  Exact k = 1/R = %.1f 1/m\n', ...
        R*1e6, 1/R);
fprintf('\n     N    dx [um]  cells/R   max err (band)   rel err    order\n');
fprintf('   ---------------------------------------------------------------\n');

for q = 1:numel(Ns)
    N = Ns(q);  dx = L/N;  dy = dx;
    xc = (dx/2:dx:L-dx/2)' - L/2;
    yc = xc;
    [X,Y] = ndgrid(xc,yc);
    phi = R - sqrt(X.^2 + Y.^2);         % >0 inside, exact signed distance

    k = curvature(phi, dx, dy);

    % Evaluate only in a narrow band around the interface -- that is the only
    % place curvature is used, and far from the interface the formula is
    % dominated by the 1/|grad phi|^3 denominator at the centre singularity.
    band = abs(phi) < 2*dx;
    e = max(abs(k(band) - 1/R));
    err(q) = e;

    if q == 1
        fprintf('   %4d   %7.3f   %6.1f   %12.4e   %8.2e     ---\n', ...
                N, dx*1e6, R/dx, e, e*R);
    else
        p = log2(err(q-1)/err(q));
        fprintf('   %4d   %7.3f   %6.1f   %12.4e   %8.2e    %5.2f\n', ...
                N, dx*1e6, R/dx, e, e*R, p);
    end
end
fprintf('\n   Expect ~2nd order while the interface is well resolved.\n');

%% ---------------------------------------------- noise sensitivity --------
fprintf('\n  Noise sensitivity (Exercise 10.1): curvature is a 2nd derivative,\n');
fprintf('  so an error of amplitude eps in phi becomes ~eps/dx^2 in k.\n\n');
N = 128;  dx = L/N;
xc = (dx/2:dx:L-dx/2)' - L/2;  [X,Y] = ndgrid(xc,xc);
phi0 = R - sqrt(X.^2+Y.^2);
band = abs(phi0) < 2*dx;

fprintf('    noise amp   max|k - 1/R|   relative to 1/R   predicted amp/dx^2\n');
fprintf('   --------------------------------------------------------------------\n');
rng(0);
for amp = [0 1e-9 1e-8 1e-7]
    phi = phi0 + amp*randn(size(phi0));
    k = curvature(phi, dx, dx);
    e = max(abs(k(band) - 1/R));
    fprintf('    %8.0e   %12.4e   %15.3f   %16.2e\n', amp, e, e*R, amp/dx^2);
end
fprintf(['\n   Noise of 1e-7 m (a tenth of a micron, 0.4%% of dx) already\n' ...
         '   swamps the true curvature. THIS is why the level set must be kept\n' ...
         '   a clean signed-distance function (Tutorial 11.3), and why a noisy\n' ...
         '   phi produces spurious currents.\n']);

%% -------------------------------------------- the pressure that implies --
fprintf('\n  Young-Laplace pressure jumps (water/air, sigma = %g N/m):\n', sigma);
fprintf('\n    diameter    2-D (cyl): sigma/R    3-D (sph): 2 sigma/R\n');
fprintf('   ------------------------------------------------------------\n');
for d = [1e-6 10e-6 100e-6 1e-3]
    Rk = d/2;
    fprintf('    %6.0f um   %14.1f Pa   %17.1f Pa\n', d*1e6, sigma/Rk, 2*sigma/Rk);
end
fprintf(['\n   Our 2-D simulations model an INFINITE CYLINDER, so use sigma/R.\n' ...
         '   Comparing a 2-D droplet against 2*sigma/R is a classic factor-of-2\n' ...
         '   error. Tutorial 11 validates against sigma/R for this reason.\n']);
fprintf(['\n   For scale: Tutorial 8 measured 0.67 Pa across an 800 um channel.\n' ...
         '   A 50 um droplet carries ~2900 Pa (2-D) of capillary pressure --\n' ...
         '   four orders of magnitude larger. That ratio is why interfaces snap\n' ...
         '   to equilibrium instantly, why a stuck bubble is so hard to move,\n' ...
         '   and why surface tension sets the timestep.\n\n']);

if ~QUIET
    N = 128; dx = L/N;
    xc = (dx/2:dx:L-dx/2)' - L/2;  [X,Y] = ndgrid(xc,xc);
    phi = R - sqrt(X.^2+Y.^2);
    k = curvature(phi,dx,dx);
    kb = k;  kb(abs(phi) > 2*dx) = NaN;

    figure('Color','w','Position',[80 80 940 380]);
    subplot(1,2,1);
    contourf(xc*1e6, xc*1e6, phi.'*1e6, 20,'LineStyle','none'); hold on
    contour(xc*1e6, xc*1e6, phi.', [0 0],'w','LineWidth',2);
    axis equal tight; colorbar; title('\phi [\mum]'); xlabel('x [\mum]');
    subplot(1,2,2);
    imagesc(xc*1e6, xc*1e6, kb.'); set(gca,'YDir','normal');
    axis equal tight; colorbar; title(sprintf('\\kappa near interface (exact %.0f)',1/R));
    xlabel('x [\mum]');
end
