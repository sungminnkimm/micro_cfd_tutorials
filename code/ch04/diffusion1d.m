%DIFFUSION1D  Explicit / implicit / Crank-Nicolson, and a TIME-order study.
%
%   Tutorial 4.  Run:  diffusion1d
%
%   Expected observed orders in time: FE 1, BE 1, CN 2.
%
%   ISOLATING THE TIME ERROR
%   ------------------------
%   Comparing against the analytic PDE solution mixes space and time error, and
%   once the time error drops below the (fixed) space error the measured order
%   collapses to zero.  Instead we compare against the exact solution of the
%   SEMI-DISCRETE system dc/dt = D*L*c, which is available in closed form here:
%   the discrete cosine is an exact eigenvector of the Neumann Laplacian.
%
%     c_j(0) = cos(k x_j),  x_j = (j-1/2)h,  k = pi/L
%     L c    = lambda c  with  lambda = -4 sin^2(k h/2) / h^2
%     => exact semi-discrete solution  c(t) = c(0) * exp(D*lambda*t)
%
%   The Neumann ghost relation c_0 = c_1 is satisfied EXACTLY by this mode
%   (x_0 = -h/2 and x_1 = +h/2 give the same cosine), so no boundary error
%   contaminates the measurement.  What remains is pure time-integration error.
%
%   See also STABILITY_DEMO, DIFFUSION2D, DIFF_MATRIX.

if ~exist('QUIET','var'), QUIET = false; end

L = 100e-6;  D = 5e-10;

%% ---------------------------------------------------------------- setup ---
N  = 100;  h = L/N;
x  = (h/2 : h : L-h/2)';
k  = pi/L;
c0 = cos(k*x);

[~, Lap] = diff_matrix(N, h, 'neumann0');
lambda   = -4*sin(k*h/2)^2 / h^2;              % exact discrete eigenvalue
fprintf('\n  discrete eigenvalue  %.6e   continuum -k^2 = %.6e\n', lambda, -k^2);
fprintf('  (they agree to O(h^2) -- that is the SPACE error, held fixed here)\n');

tEnd  = 0.02 * L^2/D;
exact = @(t) c0 * exp(D*lambda*t);             % exact semi-discrete solution

I = speye(N);

%% ------------------------------------------------- time-order study -------
% Each scheme gets its OWN starting dt.  Forward Euler cannot be run at the
% implicit schemes' timestep -- it is unstable there, and forcing it produces
% a column of NaN rather than a convergence rate.  Its dt0 sits just inside
% r = D*dt/h^2 <= 1/2; the implicit schemes start ~200x larger, which is the
% entire practical point of being unconditionally stable.
dt_stab  = 0.5*h^2/D;
dt0_all  = [0.8*dt_stab, 2e-3*L^2/D, 2e-3*L^2/D];   % FE, BE, CN
schemes  = {'FE','BE','CN'};
nlev     = 6;
err      = zeros(nlev, 3);
dtused   = zeros(nlev, 3);

for s = 1:3
    for m = 1:nlev
        dt  = dt0_all(s) / 2^(m-1);
        nst = max(1, round(tEnd/dt));
        dt  = tEnd/nst;                        % land exactly on tEnd
        dtused(m,s) = dt;

        switch schemes{s}
            case 'FE', A = [];              B = I + dt*D*Lap;
            case 'BE', A = I - dt*D*Lap;    B = I;
            case 'CN', A = I - dt*D/2*Lap;  B = I + dt*D/2*Lap;
        end
        if ~isempty(A), dA = decomposition(A,'lu'); end

        c = c0;
        for n = 1:nst
            c = B*c;
            if ~isempty(A), c = dA\c; end
        end
        err(m,s) = norm(c - exact(tEnd), inf);
    end
end

p = log2(err(1:end-1,:)./err(2:end,:));
fprintf('\n  TIME convergence (grid fixed at N = %d, t_end = %.3g s)\n', N, tEnd);
fprintf('  FE starts at dt = %.2e s (its stability limit is %.2e s)\n', ...
        dtused(1,1), dt_stab);
fprintf('  BE/CN start at dt = %.2e s -- %.0fx larger, and still stable.\n\n', ...
        dtused(1,2), dtused(1,2)/dtused(1,1));
fprintf('    level        FE err  ord |      BE err  ord |      CN err  ord\n');
fprintf('  ----------------------------------------------------------------\n');
for m = 1:nlev
    fprintf('   dt0/%-3d', 2^(m-1));
    for s = 1:3
        fprintf('  %10.3e', err(m,s));
        if m>1, fprintf(' %4.2f |', p(m-1,s)); else, fprintf('  --- |'); end
    end
    fprintf('\n');
end
fprintf('\n  Observed: FE %.2f | BE %.2f | CN %.2f   (expected 1, 1, 2)\n\n', ...
        p(end,1), p(end,2), p(end,3));

%% ------------------------------------------- physical demo: two streams ---
% Exercise 4.5's mapping: two streams meeting in a Y-junction, viewed in the
% cross-stream direction with downstream distance playing the role of time.
U     = 2e-3;                          % m/s
w     = 50e-6;                         % channel width
Nw    = 200; hw = w/Nw;
xw    = (hw/2:hw:w-hw/2)';
cw    = double(xw < w/2);              % half the channel seeded with dye
[~,Lw] = diff_matrix(Nw,hw,'neumann0');

% Crank-Nicolson is unconditionally stable, so choose dt for ACCURACY, not
% stability.  The explicit limit here would be 0.4*hw^2/D = 5e-5 s, needing
% ~100k steps; 0.5% of the diffusion time is plenty accurate and needs ~200.
dt  = 0.005 * w^2/D;
Acn = decomposition(speye(Nw) - dt*D/2*Lw, 'lu');
Bcn = speye(Nw) + dt*D/2*Lw;

c_ref  = cw;                            % unmixed reference for the index
nmax   = 2000;
hist_t = zeros(1,nmax);  hist_m = zeros(1,nmax);
hist_m(1) = mixing_index(cw, c_ref);    % = 1 by construction
n = 1;  t = 0;
while hist_m(n) > 0.05 && n < nmax
    cw = Acn\(Bcn*cw);
    t  = t + dt;
    n  = n + 1;
    hist_t(n) = t;
    hist_m(n) = mixing_index(cw, c_ref);
end
hist_t = hist_t(1:n);  hist_m = hist_m(1:n);

Pe = U*w/D;
Lmix = U*t;
fprintf('  Y-junction, w = %g um, U = %g mm/s, Pe = %.0f\n', w*1e6, U*1e3, Pe);
fprintf('    time to M < 0.05         : %.2f s\n', t);
fprintf('    => channel length needed : %.2f mm\n', Lmix*1e3);
fprintf('    Tutorial 1 rule  Pe*w    : %.2f mm  (overestimates by %.1fx)\n', ...
        Pe*w*1e3, Pe*w/Lmix);
fprintf('    calibrated rule  %.2f*Pe*w\n', Lmix/(Pe*w));
fprintf(['\n    The Pe*w rule is an ORDER-OF-MAGNITUDE estimate, not an equality.\n' ...
         '    Here the true coefficient is ~0.3, because (a) the diffusion\n' ...
         '    distance is w/2 not w, and (b) mixing is set by the decay of the\n' ...
         '    slowest Fourier mode, exp(-D pi^2 t/w^2), which carries a 1/pi^2.\n' ...
         '    Use Pe*w to decide whether a design is plausible; use a solve like\n' ...
         '    this one to size the actual channel.\n\n']);

if ~QUIET
    figure('Color','w','Position',[100 100 940 380]);
    subplot(1,2,1);
    loglog(dtused(:,1), err(:,1),'o-','LineWidth',1.5); hold on
    loglog(dtused(:,2), err(:,2),'s-','LineWidth',1.5);
    loglog(dtused(:,3), err(:,3),'^-','LineWidth',1.5);
    loglog(dtused(:,2), err(1,2)*(dtused(:,2)/dtused(1,2)),   'k--');
    loglog(dtused(:,3), err(1,3)*(dtused(:,3)/dtused(1,3)).^2,'k:');
    grid on; xlabel('\Deltat [s]'); ylabel('max-norm error');
    legend({'forward Euler','backward Euler','Crank-Nicolson','O(\Deltat)','O(\Deltat^2)'}, ...
           'Location','southeast');
    title('time-order convergence');

    subplot(1,2,2);
    plot(hist_t*U*1e3, hist_m,'LineWidth',1.6); grid on
    xlabel('downstream distance [mm]'); ylabel('mixing index (1 = unmixed)');
    title(sprintf('Y-junction mixing, Pe = %.0f', Pe));
end
