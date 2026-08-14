%CELL_PECLET  The wiggle threshold: central vs upwind on a boundary layer.
%
%   Tutorial 5.  Run:  cell_peclet
%
%   Steady 1-D advection-diffusion   u dc/dx = D d2c/dx2   on [0,L] with
%   c(0) = 0, c(L) = 1 has the exact solution
%
%       c(x) = (1 - exp(Pe*x/L)) / (1 - exp(Pe)),     Pe = u*L/D
%
%   a flat region plus an exponential boundary layer of thickness L/Pe at the
%   outflow end.  Central differencing is oscillation-free only while the CELL
%   Peclet number Pe_cell = u*dx/D <= 2 -- i.e. only while the grid can resolve
%   that layer.  Above it you get wiggles: not instability (nothing blows up),
%   just a scheme being asked to represent a layer thinner than a cell.
%
%   See also ADVECT1D, NUM_DIFFUSION.

if ~exist('QUIET','var'), QUIET = false; end

L = 1;  u = 1;  N = 50;
dx = L/N;
xf = (0:dx:L)';                       % node-centred here: BCs sit ON nodes
xi = xf(2:end-1);                     % interior unknowns
ni = numel(xi);

PeCells = [1 2 4 20];
CEN = zeros(ni, numel(PeCells));
UPW = zeros(ni, numel(PeCells));
EXA = zeros(ni, numel(PeCells));

fprintf('\n  Steady advection-diffusion, N = %d cells\n\n', N);
fprintf('   Pe_cell   Pe(global)   central: min c   wiggles?   upwind: layer err\n');
fprintf('  ---------------------------------------------------------------------\n');

for q = 1:numel(PeCells)
    Pec = PeCells(q);
    D   = u*dx/Pec;                   % pick D to hit the target cell Peclet
    Peg = u*L/D;

    e = ones(ni,1);

    % Discretize  u dc/dx - D d2c/dx2 = 0  and collect coefficients.
    % ---- central differencing
    ac_lo = -u/(2*dx) - D/dx^2;
    ac_di =             2*D/dx^2;
    ac_hi =  u/(2*dx) - D/dx^2;
    Ac = spdiags([ac_lo*e, ac_di*e, ac_hi*e], [-1 0 1], ni, ni);

    % ---- first-order upwind (u > 0, so look left)
    au_lo = -u/dx - D/dx^2;
    au_di =  u/dx + 2*D/dx^2;
    au_hi =       - D/dx^2;
    Au = spdiags([au_lo*e, au_di*e, au_hi*e], [-1 0 1], ni, ni);

    % c(0) = 0 contributes nothing;  c(L) = 1 moves to the RHS of the last row.
    bc = zeros(ni,1);  bc(end) = -ac_hi*1;
    bu = zeros(ni,1);  bu(end) = -au_hi*1;

    cc = Ac\bc;   cu = Au\bu;

    % Exact solution, written to avoid overflow.  The textbook form
    %     (1 - exp(Pe*x/L)) / (1 - exp(Pe))
    % overflows to Inf/Inf = NaN once Pe > ~709.  Multiplying top and bottom
    % by exp(-Pe) puts every exponent <= 0 and is exact for all Pe.
    ce = (exp(-Peg) - exp(Peg*(xi/L - 1))) ./ (exp(-Peg) - 1);

    CEN(:,q) = cc;  UPW(:,q) = cu;  EXA(:,q) = ce;

    wig = min(cc) < -1e-8 || max(cc) > 1+1e-8;
    if wig, wtxt = 'YES'; else, wtxt = 'no '; end
    fprintf('  %7.1f   %10.1f   %13.4f   %8s   %14.3e\n', ...
            Pec, Peg, min(cc), wtxt, norm(cu-ce,inf));
end

fprintf(['\n  Central differencing has NO numerical diffusion but wiggles once\n' ...
         '  Pe_cell > 2.  Upwind never wiggles but smears the layer badly.\n' ...
         '  At the Pe ~ 200 of a real microchannel neither is acceptable,\n' ...
         '  which is the entire motivation for flux limiters.\n\n']);

if ~QUIET
    figure('Color','w','Position',[100 100 940 640]);
    for q = 1:numel(PeCells)
        subplot(2,2,q);
        plot(xi, EXA(:,q),'k-','LineWidth',2); hold on
        plot(xi, CEN(:,q),'o-','MarkerSize',3,'LineWidth',1.1);
        plot(xi, UPW(:,q),'s-','MarkerSize',3,'LineWidth',1.1);
        grid on; xlabel('x'); ylabel('c');
        title(sprintf('Pe_{cell} = %g', PeCells(q)));
        if q==1, legend({'exact','central','upwind'},'Location','northwest'); end
    end
    sgtitle('The wiggle threshold is Pe_{cell} = 2');
end
