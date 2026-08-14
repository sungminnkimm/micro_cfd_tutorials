%STABILITY_DEMO  Watch the explicit diffusion scheme cross its stability limit.
%
%   Tutorial 4.  Run:  stability_demo
%
%   The limit r = D*dt/dx^2 <= 1/2 is SHARP.  There is no "slightly unstable";
%   there is only "growing slowly enough that you have not noticed yet".
%
%   See also DIFFUSION1D, DIFFUSION2D.

if ~exist('QUIET','var'), QUIET = false; end

L = 100e-6;  D = 5e-10;  N = 101;
dx = L/(N-1);  x = linspace(0,L,N)';
c0 = double(abs(x - L/2) < L/10);        % a top-hat of dye

rvals  = [0.40 0.51 0.60];
nsteps = 400;
C      = zeros(N, numel(rvals));

fprintf('\n  Explicit forward Euler, %d steps, stability limit r <= 0.5\n\n', nsteps);
for q = 1:numel(rvals)
    r  = rvals(q);
    c  = c0;
    for n = 1:nsteps
        c(2:end-1) = c(2:end-1) + r*(c(1:end-2) - 2*c(2:end-1) + c(3:end));
        c(1) = c(2);  c(end) = c(end-1);      % zero-flux walls
        if ~all(isfinite(c)), break; end
    end
    C(:,q) = c;

    % A sawtooth mode has large cell-to-cell curvature: measure it directly.
    % diff(c,2) is the undivided second difference -- big for (+1,-1,+1,...).
    if all(isfinite(c)), saw = max(abs(diff(c,2))); else, saw = Inf; end
    if r <= 0.5, tag = 'STABLE'; else, tag = 'UNSTABLE'; end
    fprintf('  r = %.2f : max|c| = %10.3e   sawtooth = %10.3e   %s\n', ...
            r, max(abs(c)), saw, tag);
end

fprintf(['\n  Note r = 0.51: only 1%% over the limit, still finite after %d\n' ...
         '  steps, but the sawtooth indicator has already grown. Run it for\n' ...
         '  5000 steps (Exercise 4.1) and it dies.\n\n'], nsteps);

if ~QUIET
    figure('Color','w','Position',[100 100 900 340]);
    for q = 1:numel(rvals)
        subplot(1,3,q);
        plot(x*1e6, C(:,q),'-','LineWidth',1.3); grid on
        xlabel('x [\mum]'); ylabel('c');
        title(sprintf('r = %.2f', rvals(q)));
        if all(isfinite(C(:,q))), ylim([-0.2 1.1]); end
    end
end
