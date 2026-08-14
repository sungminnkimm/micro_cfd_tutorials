%CHECKERBOARD_DEMO  Why colocated grids fail and staggered grids do not.
%
%   Tutorial 6.  Run:  checkerboard_demo
%
%   The checkerboard pressure field  p(i,j) = (-1)^(i+j)  is INVISIBLE to a
%   central-difference gradient on a colocated grid: dp/dx = (p_{i+1}-p_{i-1})/2dx
%   is identically zero because p_{i+1} and p_{i-1} are always equal.
%
%   A pressure field that exerts no force sits in the null space of the discrete
%   gradient.  The solver can add any multiple of it and nothing objects, so it
%   grows from round-off until the pressure field is garbage -- while the
%   velocity field still looks fine, which is what makes it so hard to catch.
%
%   See also MAC_GRID, POISEUILLE.

if ~exist('QUIET','var'), QUIET = false; end

nx = 16; ny = 16; dx = 1/nx; dy = 1/ny;
[I,J] = ndgrid(1:nx, 1:ny);
p = (-1).^(I+J);                      % the checkerboard mode

%% ---- colocated: central difference, needs interior nodes only ------------
dpdx_col = (p(3:end,2:end-1) - p(1:end-2,2:end-1)) / (2*dx);
dpdy_col = (p(2:end-1,3:end) - p(2:end-1,1:end-2)) / (2*dy);

%% ---- staggered: compact two-point difference AT THE FACE -----------------
dpdx_stag = (p(2:end,:) - p(1:end-1,:)) / dx;     % lives on u-faces
dpdy_stag = (p(:,2:end) - p(:,1:end-1)) / dy;     % lives on v-faces

fprintf('\n  Checkerboard pressure, nx = ny = %d, dx = %.4f\n', nx, dx);
fprintf('  amplitude of p                     : %.4f\n', max(abs(p(:))));
fprintf('\n  COLOCATED (central, p_{i+1}-p_{i-1})/2dx:\n');
fprintf('    max|dp/dx| = %.3e\n', max(abs(dpdx_col(:))));
fprintf('    max|dp/dy| = %.3e   <-- the mode is INVISIBLE\n', max(abs(dpdy_col(:))));
fprintf('\n  STAGGERED (compact, p_i - p_{i-1})/dx:\n');
fprintf('    max|dp/dx| = %.3e\n', max(abs(dpdx_stag(:))));
fprintf('    max|dp/dy| = %.3e   <-- the mode is SEEN, loudly\n', max(abs(dpdy_stag(:))));
fprintf('\n  ratio staggered/colocated = %.3e  (i.e. Inf up to round-off)\n', ...
        max(abs(dpdx_stag(:))) / max(max(abs(dpdx_col(:))), realmin));

%% ---- the same problem for velocity divergence ----------------------------
% A checkerboard VELOCITY field reports zero divergence on a colocated grid,
% so a physically impossible flow passes the incompressibility check.
u = (-1).^(I+J);  v = zeros(nx,ny);
div_col = (u(3:end,2:end-1) - u(1:end-2,2:end-1))/(2*dx);
fprintf('\n  Checkerboard VELOCITY, colocated divergence: max = %.3e\n', ...
        max(abs(div_col(:))));
fprintf('  A wildly oscillating velocity field passes as divergence-free.\n');

%% ---- the null space, counted --------------------------------------------
% Build the 1-D colocated central gradient with periodic wrap and count how
% many singular values are (numerically) zero.
e = ones(nx,1);
Gc = spdiags([-e e],[-1 1],nx,nx)/(2*dx);  Gc(1,nx) = -1/(2*dx); Gc(nx,1) = 1/(2*dx);
Gs = spdiags([-e e],[-1 0],nx,nx)/dx;      Gs(1,nx) = -1/dx;      % staggered
sc = svd(full(Gc));  ss = svd(full(Gs));
fprintf('\n  1-D periodic gradient operators, n = %d:\n', nx);
fprintf('    colocated: %d zero singular values (null space)\n', sum(sc < 1e-10*max(sc)));
fprintf('    staggered: %d zero singular values\n', sum(ss < 1e-10*max(ss)));
fprintf(['\n  The colocated operator has TWO null modes (the constant, which is\n' ...
         '  physical and harmless, plus the checkerboard, which is not). The\n' ...
         '  staggered operator has only the constant. That single extra mode is\n' ...
         '  the whole disease.\n\n']);

if ~QUIET
    figure('Color','w','Position',[100 100 900 340]);
    subplot(1,3,1); imagesc(p.'); axis equal tight; colorbar;
    set(gca,'YDir','normal'); title('checkerboard p');
    subplot(1,3,2); imagesc(dpdx_col.'); axis equal tight; colorbar;
    set(gca,'YDir','normal'); title('colocated dp/dx  (all zero)');
    subplot(1,3,3); imagesc(dpdx_stag.'); axis equal tight; colorbar;
    set(gca,'YDir','normal'); title('staggered dp/dx  (large)');
end
