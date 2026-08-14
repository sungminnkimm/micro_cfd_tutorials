%ADVECT1D  Compare advection schemes on a periodic domain.
%
%   Tutorial 5.  Run:  advect1d
%
%   A top hat and a Gaussian are advected once around a periodic domain, so the
%   exact answer at the end is the initial condition.  Anything that differs is
%   pure numerical error, and you can see exactly what kind:
%     central  -> unconditionally unstable, blows up
%     upwind   -> stable, but smears (numerical diffusion)
%     vanleer  -> stable AND sharp
%     superbee -> sharpest, slightly over-steepens a smooth Gaussian
%
%   See also LIMITER, CELL_PECLET, NUM_DIFFUSION.

if ~exist('QUIET','var'), QUIET = false; end

L  = 1;  u = 1;  N = 200;
dx = L/N;
x  = (dx/2 : dx : L-dx/2)';
C  = 0.5;                                   % Courant number
dt = C*dx/abs(u);
nsteps = round(L/abs(u)/dt);                % exactly one lap
dt = (L/abs(u))/nsteps;  C = u*dt/dx;

tophat = double(abs(x-0.3) < 0.1);
gauss  = exp(-((x-0.7).^2)/(2*0.04^2));
c0     = tophat + gauss;

schemes = {'central','upwind','minmod','vanleer','superbee','lw'};
CC = zeros(N, numel(schemes));

fprintf('\n  One lap of a periodic domain, N = %d, C = %.2f, %d steps\n\n', ...
        N, C, nsteps);
fprintf('  scheme       peak(hat)  peak(gauss)    L1 error      min(c)\n');
fprintf('  --------------------------------------------------------------\n');

for s = 1:numel(schemes)
    c = c0;
    for n = 1:nsteps
        if strcmp(schemes{s},'central')
            % Unconditionally unstable -- included so you can watch it fail.
            c = c - C/2*(circshift(c,-1) - circshift(c,1));
        else
            c = advect_tvd(c, u, dt, dx, schemes{s});
        end
        if ~all(isfinite(c)), break; end
    end
    CC(:,s) = c;

    if all(isfinite(c))
        pk_h = max(c(x<0.5));  pk_g = max(c(x>=0.5));
        l1   = sum(abs(c-c0))*dx;
        und  = min(c);
        fprintf('  %-10s  %10.4g   %10.4g   %10.3e  %+10.2e\n', ...
                schemes{s}, pk_h, pk_g, l1, und);
    else
        fprintf('  %-10s  %10s   %10s   %10s  %10s\n', ...
                schemes{s}, 'DIVERGED','-','-','-');
    end
end

fprintf(['\n  Exact answer: peak(hat) = 1, peak(gauss) = 1, L1 error = 0,\n' ...
         '  min(c) = 0 (any negative value is an unphysical undershoot).\n\n' ...
         '  Read the table as a trade-off, not a ranking:\n' ...
         '    central   diverges -- unconditionally unstable, as derived.\n' ...
         '    upwind    never undershoots but loses the most amplitude.\n' ...
         '    lw        keeps amplitude but undershoots badly (the wiggles),\n' ...
         '              and overshoots the top hat to 1.23.\n' ...
         '    TVD       most of lw''s sharpness with NO undershoot. That is\n' ...
         '              the entire reason flux limiters exist.\n' ...
         '  Note superbee beats vanleer on the top hat but slightly squares\n' ...
         '  off the smooth Gaussian -- sharpening is not free.\n\n']);

if ~QUIET
    figure('Color','w','Position',[100 100 940 560]);
    for s = 1:numel(schemes)
        subplot(2,3,s);
        plot(x, c0,'k--','LineWidth',1); hold on
        if all(isfinite(CC(:,s)))
            plot(x, CC(:,s),'-','LineWidth',1.5);
            ylim([-0.3 1.3]);
        else
            text(0.5,0.5,'DIVERGED','HorizontalAlignment','center', ...
                 'Color','r','FontWeight','bold');
        end
        grid on; title(schemes{s}); xlabel('x');
        if s==1, legend({'exact','computed'},'Location','north'); end
    end
    sgtitle(sprintf('One lap, C = %.2f', C));
end
