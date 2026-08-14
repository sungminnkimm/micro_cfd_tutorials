%FD_STENCILS  Visual sanity check of the derivative matrices.
%
%   Tutorial 2.  Run:  fd_stencils
%
%   Two things every CFD programmer should do reflexively before trusting a new
%   operator:  (1) print it small and read the numbers, (2) apply it to a
%   function whose derivative you know by hand.
%
%   See also DIFF_MATRIX, CONVERGENCE_STUDY.

if ~exist('QUIET','var'), QUIET = false; end

%% (1) Read the matrix
N = 8; h = 1/N;
fprintf('\n--- D2, dirichlet0, N=%d, h=%.4g ---\n', N, h);
[~, D2d] = diff_matrix(N,h,'dirichlet0');
disp(full(D2d)*h^2);   % scaled by h^2 so the integers are readable
fprintf('(shown scaled by h^2.  Note the -3 in the corners: ghost u_0 = -u_1)\n');

fprintf('\n--- D2, neumann0 ---\n');
[~, D2n] = diff_matrix(N,h,'neumann0');
disp(full(D2n)*h^2);
fprintf('(corners are -1: ghost u_0 = +u_1, zero flux through the wall)\n');

%% (2) Apply to a known function
% Cell-centred grid on [0,1]:  x_j = (j-1/2)h
Nf = 200; hf = 1/Nf;
x  = (hf/2 : hf : 1-hf/2)';

% u = sin(pi x) vanishes at BOTH walls, so it is compatible with dirichlet0.
u    = sin(pi*x);
d2ex = -pi^2*sin(pi*x);
[~,D2] = diff_matrix(Nf,hf,'dirichlet0');
err = norm(D2*u - d2ex, inf)/norm(d2ex, inf);
fprintf('\nDirichlet check  u=sin(pi x):  relative max error = %.3e\n', err);

% u = cos(pi x) has zero slope at BOTH walls -> compatible with neumann0.
u    = cos(pi*x);
d2ex = -pi^2*cos(pi*x);
[~,D2] = diff_matrix(Nf,hf,'neumann0');
err = norm(D2*u - d2ex, inf)/norm(d2ex, inf);
fprintf('Neumann   check  u=cos(pi x):  relative max error = %.3e\n', err);

% Mismatch demo: apply the Dirichlet operator to a function that does NOT
% vanish at the walls.  The interior is fine; the boundary rows are garbage.
u    = cos(pi*x);
d2ex = -pi^2*cos(pi*x);
[~,D2] = diff_matrix(Nf,hf,'dirichlet0');
e = abs(D2*u - d2ex);
fprintf(['\nMISMATCHED bc (Dirichlet operator, non-vanishing function):\n' ...
         '   interior max error = %.3e   <-- fine\n' ...
         '   boundary   error   = %.3e   <-- 100%% wrong\n'], ...
         max(e(3:end-2)), max(e([1 end])));
fprintf(['This is what a wrong boundary condition looks like: the field is\n' ...
         'correct everywhere you look first, and wrong only in the two cells\n' ...
         'you never plot.\n\n']);

if ~QUIET
    figure('Color','w','Position',[100 100 900 380]);
    subplot(1,3,1); spy(D2d);  title('D2 sparsity (tridiagonal)');
    subplot(1,3,2);
    [~,D2p] = diff_matrix(N,h,'periodic'); spy(D2p);
    title('periodic: note the corners');
    subplot(1,3,3);
    semilogy(x, max(e,1e-18),'LineWidth',1.4); grid on
    xlabel('x'); ylabel('|error|'); title('mismatched BC: error is at the walls');
end
