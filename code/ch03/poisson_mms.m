%POISSON_MMS  Verify poisson2d by the method of manufactured solutions.
%
%   Tutorial 3.  Run:  poisson_mms
%
%   Pick the answer, differentiate it to get the source, feed the source to the
%   solver, compare.  Works for any operator and any BC combination -- this is
%   the professional standard for code verification (Tutorial 13).
%
%   See also POISSON2D, POISSON_SOLVE.

if ~exist('QUIET','var'), QUIET = false; end

phi_ex = @(x,y) sin(pi*x).*cos(2*pi*y);
src    = @(x,y) -(pi^2 + 4*pi^2)*sin(pi*x).*cos(2*pi*y);   % = del^2(phi_ex)

N   = [16 32 64 128 256];
err = zeros(size(N));

for k = 1:numel(N)
    % Deliberately NON-SQUARE: catches kron-ordering and ndgrid/meshgrid bugs
    % that are invisible when nx == ny.
    nx = N(k);  ny = round(N(k)*3/4);
    dx = 1/nx;  dy = 1/ny;
    xc = (dx/2 : dx : 1-dx/2)';        % cell centres
    yc = (dy/2 : dy : 1-dy/2)';
    [X,Y] = ndgrid(xc,yc);             % ndgrid: index order (x,y)

    S = poisson2d(nx,ny,dx,dy,{'D','D','D','D'});

    % Dirichlet values evaluated ON the walls (x = 0,1 and y = 0,1)
    bcVals = { phi_ex(0,yc),  phi_ex(1,yc), ...
               phi_ex(xc,0),  phi_ex(xc,1) };

    phi    = poisson_solve(S, src(X,Y), bcVals);
    err(k) = norm(phi(:) - reshape(phi_ex(X,Y),[],1), inf);
end

p = log2(err(1:end-1)./err(2:end));
fprintf('\n    nx    max error   order\n');
fprintf('  ------------------------------\n');
for k = 1:numel(N)
    if k==1, fprintf('  %4d   %9.3e     ---\n', N(k), err(k));
    else,    fprintf('  %4d   %9.3e    %5.2f\n', N(k), err(k), p(k-1)); end
end
fprintf('\n  Observed order: %.2f  (expected 2.00)\n\n', p(end));

%% Neumann check: does the compatibility condition hold?
% All-Neumann with phi = cos(pi x)cos(pi y): dphi/dn = 0 on every wall, and
% int(f)dV = 0 automatically.  The pinned solve recovers phi up to a constant,
% so compare AFTER removing the mean.
nx = 64; ny = 48; dx = 1/nx; dy = 1/ny;
xc = (dx/2:dx:1-dx/2)'; yc = (dy/2:dy:1-dy/2)';
[X,Y] = ndgrid(xc,yc);
pe = cos(pi*X).*cos(pi*Y);
fe = -2*pi^2*pe;

Sn  = poisson2d(nx,ny,dx,dy,{'N','N','N','N'});
phi = poisson_solve(Sn, fe, {0,0,0,0});

phi = phi - mean(phi(:));            % remove the arbitrary constant
pe  = pe  - mean(pe(:));
fprintf('All-Neumann (pinned), mean-removed max error: %.3e\n', ...
        norm(phi(:)-pe(:),inf));
fprintf('Compatibility  int(f)dV = %.3e  (must be ~0 for a solution to exist)\n\n', ...
        sum(fe(:))*dx*dy);

if ~QUIET
    figure('Color','w','Position',[100 100 900 350]);
    subplot(1,3,1); contourf(xc,yc,pe.',20,'LineStyle','none');
    axis equal tight; colorbar; title('exact');
    subplot(1,3,2); contourf(xc,yc,phi.',20,'LineStyle','none');
    axis equal tight; colorbar; title('computed');
    subplot(1,3,3); contourf(xc,yc,(phi-pe).',20,'LineStyle','none');
    axis equal tight; colorbar; title('error');
end
