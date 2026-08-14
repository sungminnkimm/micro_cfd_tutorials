%MENISCUS  Arc shapes, capillary entry pressure, and a burst valve.
%
%   Tutorial 10.  Run:  meniscus
%
%   At Bo << 1 gravity is negligible, so a static meniscus has CONSTANT
%   curvature -- a circular arc meeting each wall at the contact angle.  For a
%   2-D channel of width w:
%
%       R  = w / (2 cos theta)
%       dp = sigma / R = 2 sigma cos(theta) / w      <- capillary entry pressure
%
%   That dp is the pressure needed to push an interface through a constriction
%   of width w, and it governs spontaneous wicking, burst valves, and why a
%   narrow neck blocks a non-wetting phase.
%
%   See also CAPILLARY_RISE, YOUNG_LAPLACE.

if ~exist('QUIET','var'), QUIET = false; end

sigma = 0.072;

%% ------------------------------------------- entry pressure vs geometry --
fprintf('\n  Capillary entry pressure  dp = 2 sigma cos(theta)/w   [mbar]\n');
fprintf('  (positive = pulls itself in;  negative = must be forced)\n\n');
ws  = [10e-6 20e-6 50e-6 100e-6 200e-6];
ths = [0 45 90 105 135];

fprintf('   w [um] |');  fprintf('  th=%3d', ths);  fprintf('\n');
fprintf('  --------+');  fprintf('%s', repmat('-',1,8*numel(ths)));  fprintf('\n');
for i = 1:numel(ws)
    fprintf('   %5.0f  |', ws(i)*1e6);
    for j = 1:numel(ths)
        dp = 2*sigma*cosd(ths(j))/ws(i);
        fprintf(' %7.1f', dp/100);       % Pa -> mbar
    end
    fprintf('\n');
end
fprintf(['\n  theta = 90 gives exactly zero: a neutral surface neither pulls the\n' ...
         '  liquid in nor pushes it out, at any width.\n']);

%% ------------------------------------------------------ worked example ---
w = 20e-6;  th = 105;
dp = 2*sigma*cosd(th)/w;
fprintf('\n  Exercise 10.4: water into a %g um NATIVE PDMS constriction (theta=%d)\n', ...
        w*1e6, th);
fprintf('    entry pressure = %.0f Pa = %.0f mbar  (negative: it resists)\n', dp, dp/100);
fprintf('    so you must apply +%.0f mbar to push water in.\n', -dp/100);
fprintf('    A syringe pump easily reaches this; the limit is usually the\n');
fprintf('    bond strength of the chip, not the pump.\n');

%% ---------------------------------------------------------- burst valve --
fprintf('\n  Exercise 10.5: burst valve by sudden expansion\n');
fprintf('  A hydrophilic channel (theta = 30) wicks spontaneously. At a sudden\n');
fprintf('  expansion the interface must bulge to a much larger radius, which\n');
fprintf('  costs pressure. Flow halts until the burst pressure is applied.\n\n');
th_h = 30;  w1 = 50e-6;
dp_in = 2*sigma*cosd(th_h)/w1;
fprintf('    narrow section w1 = %g um: wicking pressure = %+.0f mbar\n', ...
        w1*1e6, dp_in/100);
fprintf('\n    expansion ratio   w2 [um]   holding dp [mbar]   burst target 50 mbar?\n');
fprintf('   ---------------------------------------------------------------------\n');
for ratio = [2 4 8 16 32]
    w2 = w1*ratio;
    % At the expansion the interface must locally reach the maximum curvature
    % the corner allows; the classic estimate uses an effective angle of
    % theta + the expansion half-angle, which for a SHARP (90 deg) step gives
    % an effective contact angle of theta + 90.
    th_eff = min(th_h + 90, 180);
    dp_hold = -2*sigma*cosd(th_eff)/w1;     % resisting pressure at the step
    if dp_hold/100 >= 50, ok = 'yes'; else, ok = 'no'; end
    fprintf('   %10g       %7.0f   %14.0f      %s\n', ratio, w2*1e6, dp_hold/100, ok);
end
fprintf(['\n    Note the holding pressure is set by the SHARP CORNER (through the\n' ...
         '    effective angle theta+90) and by w1 -- NOT by the expansion ratio.\n' ...
         '    Making the downstream channel wider does not make the valve\n' ...
         '    stronger; making the upstream channel NARROWER does. That is the\n' ...
         '    useful design insight, and it is easy to get backwards.\n']);
fprintf('    Required w1 for a 50 mbar burst: %.1f um\n', ...
        -2*sigma*cosd(min(th_h+90,180))/(50*100)*1e6);

%% ------------------------------------------------------- arc geometry ----
if ~QUIET
    w0 = 100e-6;
    figure('Color','w','Position',[80 80 900 380]);
    subplot(1,2,1); hold on
    cols = lines(5);
    for j = 1:numel(ths)
        th = ths(j);
        if abs(cosd(th)) < 1e-6, continue; end
        R = w0/(2*cosd(th));
        yy = linspace(-w0/2, w0/2, 200);
        xx = sign(R)*(sqrt(max(R^2 - yy.^2,0)) - abs(R)*cosd(0));
        xx = sqrt(max(R^2-yy.^2,0));
        xx = sign(R)*(xx - max(xx));
        plot(xx*1e6, yy*1e6, 'LineWidth',2,'Color',cols(j,:), ...
             'DisplayName',sprintf('\\theta = %d^\\circ',th));
    end
    yline( w0/2*1e6,'k','LineWidth',2,'HandleVisibility','off');
    yline(-w0/2*1e6,'k','LineWidth',2,'HandleVisibility','off');
    axis equal; grid on; xlabel('x [\mum]'); ylabel('y [\mum]');
    legend('Location','best'); title(sprintf('Meniscus arcs, w = %g \\mum', w0*1e6));

    subplot(1,2,2);
    wv = logspace(-6,-3,200);
    for j = 1:numel(ths)
        semilogx(wv*1e6, 2*sigma*cosd(ths(j))./wv/100,'LineWidth',1.8); hold on
    end
    yline(0,'k--'); grid on;
    xlabel('channel width [\mum]'); ylabel('entry pressure [mbar]');
    legend(arrayfun(@(t)sprintf('\\theta=%d',t),ths,'UniformOutput',false));
    title('capillary entry pressure');
end
fprintf('\n');
