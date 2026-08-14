%SCALE_SWEEP  How the dimensionless groups vary with device size.
%
%   Tutorial 1.  Run:  scale_sweep
%
%   See also SCALING, SCALING_DEMO.

if ~exist('QUIET','var'), QUIET = false; end

L = logspace(-7, -2, 200);      % 100 nm .. 1 cm
U = 1e-3;                       % fixed 1 mm/s

Re = zeros(size(L)); Pe = Re; Ca = Re; Bo = Re;
for k = 1:numel(L)
    s = scaling('L',L(k),'U',U);
    Re(k)=s.Re; Pe(k)=s.Pe; Ca(k)=s.Ca; Bo(k)=s.Bo;
end

% Crossover sizes (where each group = 1).  Interpolate in log space: each group
% is a pure power of L, so log(group) is linear in log(L) and this is exact.
L_Re = exp(interp1(log(Re), log(L), 0, 'linear', NaN));
L_Pe = exp(interp1(log(Pe), log(L), 0, 'linear', NaN));
L_Bo = exp(interp1(log(Bo), log(L), 0, 'linear', NaN));

fprintf('At U = %g mm/s:\n', U*1e3);
fprintf('  Re = 1 at L = %8.3g mm  (above: inertia matters)\n', L_Re*1e3);
fprintf('  Pe = 1 at L = %8.3g um  (above: advection beats diffusion)\n', L_Pe*1e6);
fprintf('  Bo = 1 at L = %8.3g mm  (above: gravity beats surface tension)\n', L_Bo*1e3);

if ~QUIET
    figure('Color','w','Position',[100 100 760 460]);
    xr = [10 500];   % typical microfluidics band, in um
    patch([xr fliplr(xr)], [1e-10 1e-10 1e8 1e8], ...
          [0.2 0.5 0.9],'FaceAlpha',0.07,'EdgeColor','none'); hold on
    loglog(L*1e6, Re,'LineWidth',2);
    loglog(L*1e6, Pe,'LineWidth',2);
    loglog(L*1e6, Ca,'LineWidth',2);
    loglog(L*1e6, Bo,'LineWidth',2);
    yline(1,'k--','= 1','LineWidth',1.2);
    set(gca,'XScale','log','YScale','log');
    text(70, 1e6, 'typical microfluidics','Color',[0.2 0.4 0.7]);
    grid on; xlabel('channel size L  [\mum]'); ylabel('dimensionless group');
    legend({'','Re','Pe','Ca','Bo'},'Location','northwest');
    title(sprintf('Scaling at U = %g mm/s (water)', U*1e3));
    ylim([1e-10 1e8]); xlim([0.1 1e4]);
end
