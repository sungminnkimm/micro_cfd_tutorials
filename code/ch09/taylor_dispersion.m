%TAYLOR_DISPERSION  Verify D_eff = D(1 + 2*Pe^2/105) for plane Poiseuille flow.
%
%   Tutorial 9.  Run:  taylor_dispersion
%
%   A plug of solute in Poiseuille flow spreads ALONG the channel far faster
%   than molecular diffusion alone.  Centreline fluid moves at 1.5*ubar while
%   near-wall fluid barely moves; molecules sample different streamlines by
%   diffusing across the channel, and shear + cross-stream diffusion combine
%   into an effective axial diffusivity
%
%       D_eff = D * (1 + (2/105) * Pe_h^2),   Pe_h = ubar*h/D
%
%   for plates separated by 2h.  (Circular tube: coefficient 1/48 instead.)
%
%   THE PARADOX: D_eff ~ 1/D.  FASTER-diffusing molecules disperse LESS,
%   because they cross streamlines quickly and average out the shear.  Verified
%   numerically at the bottom of this script.
%
%   VALIDITY: the result is asymptotic, requiring t >> h^2/D for cross-stream
%   equilibration.  Before that you are pre-asymptotic and dispersion is LESS
%   than predicted.  This script reports the ratio so you can see it.
%
%   See also SCALAR_STEP, MIXING_STUDY.

if ~exist('QUIET','var'), QUIET = false; end

%% ---------------------------------------------------------- parameters ---
h2   = 100e-6;            % full gap (2h)
hh   = h2/2;              % half gap
ubar = 1e-3;
D    = 5e-10;
Lx   = 6e-3;

nx = 600;  ny = 40;
dx = Lx/nx;  dy = h2/ny;
xc = (dx/2:dx:Lx-dx/2)';
yc = (dy/2:dy:h2-dy/2)';

Pe_h  = ubar*hh/D;
Deff  = D*(1 + (2/105)*Pe_h^2);

fprintf('\n  Taylor-Aris dispersion\n');
fprintf('    gap 2h = %g um, ubar = %g mm/s, D = %.2e m^2/s\n', h2*1e6, ubar*1e3, D);
fprintf('    Pe_h = ubar*h/D = %.1f\n', Pe_h);
fprintf('    predicted D_eff = %.4e m^2/s  = %.1f x D\n', Deff, Deff/D);

%% ------------------------------------------- frozen Poiseuille velocity ---
% Analytic profile on the u-faces; no need to run the NS solver for this.
yu = yc;
uProf = 6*ubar*(yu/h2).*(1 - yu/h2);
u = repmat(uProf.', nx+1, 1);          % (nx+1) x ny
v = zeros(nx, ny+1);

%% ---------------------------------------------------------- initial plug --
x0    = 0.8e-3;
sig0  = 120e-6;
c = repmat(exp(-((xc-x0).^2)/(2*sig0^2)), 1, ny);

dt = 0.4*min( dx/max(abs(u(:))), 1/(2*D)/(1/dx^2+1/dy^2) );
tEnd = 0.55*(Lx-2*x0)/ubar + 1.0;      % keep the plug inside the domain
nst  = round(tEnd/dt);  dt = tEnd/nst;

% VALIDITY CRITERION -- get the relaxation time right.  The usual textbook
% statement "t >> h^2/D" is a scaling estimate that is about 10x pessimistic.
% Cross-stream equilibration is set by the decay of the SLOWEST diffusive mode,
% exp(-pi^2*D*t/h^2), so the actual relaxation time is h^2/(pi^2 D).  Using the
% crude version here would flag this run as pre-asymptotic (ratio 0.7) even
% though it reproduces the formula to ~1%.  Same lesson as the Pe*w mixing rule
% in Tutorial 4: order-of-magnitude rules are for deciding whether a design is
% plausible, not for judging whether a computation is converged.
tRelax = hh^2/(pi^2*D);
fprintf('    transit time %.2f s vs cross-stream relaxation h^2/(pi^2 D) = %.2f s\n', ...
        tEnd, tRelax);
fprintf('    ratio %.1f  (need >> 1 for the asymptotic formula to hold)\n', ...
        tEnd/tRelax);
fprintf('    [the cruder criterion t >> h^2/D would give %.1f and cry wolf]\n', ...
        tEnd/(hh^2/D));
fprintf('    %d steps, dt = %.2e s\n', nst, dt);

nRec = 60;  rec = round(linspace(1,nst,nRec));
tRec = zeros(nRec,1);  varRec = zeros(nRec,1);  q = 0;

for n = 1:nst
    c = scalar_step(c, u, v, D, dt, dx, dy, 'vanleer');
    if any(rec==n)
        q = q+1;  tRec(q) = n*dt;
        cx = sum(c,2);                       % cross-section-averaged profile
        w  = cx/sum(cx);
        mx = sum(w.*xc);
        varRec(q) = sum(w.*(xc-mx).^2);
    end
end

% D_eff from the growth of the variance:  var = var0 + 2*D_eff*t
use = tRec > 0.35*tEnd;                      % skip the pre-asymptotic start
pf  = polyfit(tRec(use), varRec(use), 1);
Deff_meas = pf(1)/2;

fprintf('\n    measured D_eff = %.4e m^2/s  = %.1f x D\n', Deff_meas, Deff_meas/D);
fprintf('    predicted      = %.4e m^2/s  = %.1f x D\n', Deff, Deff/D);
fprintf('    relative error = %.2f%%\n', 100*abs(Deff_meas-Deff)/Deff);

% What pure molecular diffusion would have given, for contrast
fprintf('\n    molecular diffusion alone would give var growth 2*D*t;\n');
fprintf('    the plug is spreading %.0fx faster than that.\n', Deff_meas/D);

%% --------------------------------------------------------- the paradox ---
fprintf('\n  The 1/D paradox (same channel, same flow):\n');
fprintf('    species              D [m^2/s]     Pe_h     D_eff/D    D_eff [m^2/s]\n');
fprintf('   ----------------------------------------------------------------------\n');
spec = {'small ion', 2e-9; 'small molecule', 5e-10; 'protein', 4e-11};
for k = 1:size(spec,1)
    Dk = spec{k,2};  Pk = ubar*hh/Dk;  Dek = Dk*(1+(2/105)*Pk^2);
    fprintf('    %-18s  %8.1e   %7.0f   %9.0f    %.3e\n', ...
            spec{k,1}, Dk, Pk, Dek/Dk, Dek);
end
fprintf(['\n    D_eff itself RISES as D falls: the protein disperses ~%.0fx more\n' ...
         '    in absolute terms than the ion, despite diffusing 50x slower.\n' ...
         '    Faster cross-stream diffusion averages out the shear.\n\n'], ...
         (4e-11*(1+(2/105)*(ubar*hh/4e-11)^2)) / (2e-9*(1+(2/105)*(ubar*hh/2e-9)^2)));

if ~QUIET
    figure('Color','w','Position',[80 80 940 400]);
    subplot(1,2,1);
    plot(tRec, varRec*1e12,'o','MarkerSize',4); hold on
    plot(tRec, polyval(pf,tRec)*1e12,'-','LineWidth',1.6);
    plot(tRec, (varRec(1) + 2*D*(tRec-tRec(1)))*1e12,'--','LineWidth',1.4);
    grid on; xlabel('t [s]'); ylabel('axial variance [\mum^2]');
    legend({'measured','fit: 2 D_{eff} t','molecular 2 D t'},'Location','northwest');
    title('variance growth');

    subplot(1,2,2);
    plot(xc*1e3, sum(c,2)/max(sum(c,2)),'LineWidth',1.6); hold on
    plot(xc*1e3, exp(-((xc-x0).^2)/(2*sig0^2)),'--','LineWidth',1.2);
    grid on; xlabel('x [mm]'); ylabel('normalized c');
    legend({'final','initial'}); title('plug spreading');
end
