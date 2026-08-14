%CONVERGENCE_STUDY  Measure the observed order of accuracy of FD stencils.
%
%   Tutorial 2.  Run:  convergence_study
%
%   The single most important diagnostic in CFD: solve a problem whose answer
%   you know, on a sequence of refined grids, and check that the error falls at
%   the theoretical rate.  Expected orders here: 1, 2, 2.
%
%   See also DIFF_MATRIX.

if ~exist('QUIET','var'), QUIET = false; end

% Test function on [0,1] with known derivatives.  Choose something with all
% derivatives nonzero so no error term accidentally vanishes and flatters us.
f   = @(x) sin(4*x) .* exp(-x);
df  = @(x) 4*cos(4*x).*exp(-x) - sin(4*x).*exp(-x);
d2f = @(x) -16*sin(4*x).*exp(-x) - 8*cos(4*x).*exp(-x) + sin(4*x).*exp(-x);

N = 2.^(4:11);                 % 16 .. 2048 intervals
h = 1./N;
[e_fwd, e_cen, e_2nd] = deal(zeros(size(N)));

for k = 1:numel(N)
    x = linspace(0,1,N(k)+1)';
    u = f(x);

    % Interior nodes only: the boundaries need one-sided stencils, and mixing
    % those in would pollute the very order we are trying to measure.
    i = 2:N(k);

    d_fwd = (u(i+1) - u(i)) / h(k);                    % O(h)
    d_cen = (u(i+1) - u(i-1)) / (2*h(k));              % O(h^2)
    d_2nd = (u(i+1) - 2*u(i) + u(i-1)) / h(k)^2;       % O(h^2)

    % Max-norm error.  (Use inf-norm; the plain 2-norm needs a sqrt(h) weight
    % to be a true function norm -- a classic way to fake a good-looking rate.)
    e_fwd(k) = norm(d_fwd - df(x(i)),  inf);
    e_cen(k) = norm(d_cen - df(x(i)),  inf);
    e_2nd(k) = norm(d_2nd - d2f(x(i)), inf);
end

% Observed order between successive grids: p = log2(E(h)/E(h/2))
p_fwd = log2(e_fwd(1:end-1)./e_fwd(2:end));
p_cen = log2(e_cen(1:end-1)./e_cen(2:end));
p_2nd = log2(e_2nd(1:end-1)./e_2nd(2:end));

fprintf('\n     N   fwd err   order |  cen err   order |  2nd err   order\n');
fprintf('  ------------------------------------------------------------\n');
for k = 1:numel(N)
    fprintf('  %4d  %9.2e', N(k), e_fwd(k));
    if k>1, fprintf(' %6.2f |', p_fwd(k-1)); else, fprintf('    --- |'); end
    fprintf(' %9.2e', e_cen(k));
    if k>1, fprintf(' %6.2f |', p_cen(k-1)); else, fprintf('    --- |'); end
    fprintf(' %9.2e', e_2nd(k));
    if k>1, fprintf(' %6.2f\n', p_2nd(k-1)); else, fprintf('    ---\n'); end
end
fprintf('\n  Asymptotic orders: fwd %.2f | cen %.2f | 2nd %.2f\n', ...
        p_fwd(end), p_cen(end), p_2nd(end));
fprintf('  Expected:          fwd 1.00 | cen 2.00 | 2nd 2.00\n\n');

if ~QUIET
    figure('Color','w','Position',[100 100 620 480]);
    loglog(h, e_fwd,'o-','LineWidth',1.6); hold on
    loglog(h, e_cen,'s-','LineWidth',1.6);
    loglog(h, e_2nd,'^-','LineWidth',1.6);
    loglog(h, h,    'k--');
    loglog(h, h.^2, 'k:');
    grid on; xlabel('h = \Deltax'); ylabel('max-norm error');
    legend({'forward d/dx','central d/dx','central d^2/dx^2','O(h)','O(h^2)'}, ...
           'Location','southeast');
    title('Grid convergence: slope = observed order');
end
