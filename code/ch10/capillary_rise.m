%CAPILLARY_RISE  How strong is capillarity? Strong enough to lift water 29 cm.
%
%   Tutorial 10.  Run:  capillary_rise
%
%       h = 2 sigma cos(theta) / (rho g r)
%
%   from a force balance: upward 2*pi*r*sigma*cos(theta) against the weight
%   rho*g*pi*r^2*h of the raised column.
%
%   See also MENISCUS, YOUNG_LAPLACE.

if ~exist('QUIET','var'), QUIET = false; end

sigma = 0.072;  rho = 998;  g = 9.81;
lc = sqrt(sigma/(rho*g));

fprintf('\n  Capillary rise of water (sigma = %g N/m)\n', sigma);
fprintf('  capillary length lc = sqrt(sigma/rho g) = %.2f mm\n\n', lc*1e3);

fprintf('   radius     theta=0     theta=30    theta=60    theta=105 (PDMS)\n');
fprintf('  ---------------------------------------------------------------------\n');
radii = [5e-6 25e-6 50e-6 250e-6 1e-3 5e-3];
thetas = [0 30 60 105];
H = zeros(numel(radii), numel(thetas));
for i = 1:numel(radii)
    fprintf('  %7.0f um', radii(i)*1e6);
    for j = 1:numel(thetas)
        H(i,j) = 2*sigma*cosd(thetas(j))/(rho*g*radii(i));
        fprintf('  %10s', fmt_h(H(i,j)));
    end
    fprintf('\n');
end

fprintf(['\n  Negative = the liquid is PUSHED DOWN. Water in native PDMS\n' ...
         '  (theta = 105 deg) will not wick into a channel at all; it must be\n' ...
         '  forced. Oxygen-plasma treatment flips the sign, which is why that\n' ...
         '  step exists in every PDMS protocol -- and why it matters that the\n' ...
         '  treatment decays over hours to days.\n']);

%% ---------------------------------------------------- the 1/r scaling ----
r = logspace(-6,-2,200);
h0 = 2*sigma./(rho*g*r);
slope = polyfit(log(r), log(h0), 1);
fprintf('\n  Scaling check: d(log h)/d(log r) = %.4f  (exactly -1)\n', slope(1));

% Where does the rise equal the capillary length?
rEq = 2*sigma/(rho*g*lc);
fprintf('  rise = capillary length at r = %.3f mm\n', rEq*1e3);
fprintf('  (i.e. h = lc when r = 2*lc -- the two lengths cross at the same\n');
fprintf('   scale, which is exactly what "lc is THE length scale" means)\n');

%% ------------------------------------------------ 100 um channel, real ---
r0 = 50e-6;
fprintf('\n  Worked example: %g um radius glass capillary, theta ~ 0\n', r0*1e6);
fprintf('    rise = %.3f m = %.1f cm\n', 2*sigma/(rho*g*r0), 2*sigma/(rho*g*r0)*100);
fprintf('    Bo = (r/lc)^2 = %.2e  -- gravity is utterly negligible ACROSS the\n', ...
        (r0/lc)^2);
fprintf('    channel; it only matters over the metre-scale rise height.\n\n');

if ~QUIET
    figure('Color','w','Position',[80 80 620 460]);
    loglog(r*1e6, h0*1e3,'LineWidth',2); hold on
    yline(lc*1e3,'k--','capillary length');
    xline(r0*1e6,'r:','50 \mum');
    grid on; xlabel('capillary radius [\mum]'); ylabel('rise height [mm]');
    title('Capillary rise, water, \theta = 0');
end

function s = fmt_h(h)
if abs(h) >= 1
    s = sprintf('%.2f m', h);
elseif abs(h) >= 1e-2
    s = sprintf('%.1f cm', h*100);
else
    s = sprintf('%.2f mm', h*1e3);
end
end
