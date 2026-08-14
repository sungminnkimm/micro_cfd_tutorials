%REGIME_MAP  Where does my droplet design sit? (No simulation required.)
%
%   Tutorial 12.  Run:  regime_map
%
%   RUN THIS BEFORE COMMITTING HOURS OF CPU.  It evaluates the dimensionless
%   groups, the Garstecki size prediction, and -- crucially -- the timestep and
%   step count your grid would demand.  It will often tell you the simulation
%   you were about to launch needs 1e7 steps and that you should change the
%   grid or the fluids first.
%
%   It also works the Tutorial 12.7 design problem, which is deliberately
%   SELF-INCONSISTENT: the naive parameter choice lands in dripping when
%   squeezing was required.  Finding that out costs nothing here.
%
%   See also T_JUNCTION, FLOW_FOCUSING, NS2PHASE.

if ~exist('QUIET','var'), QUIET = false; end

%% ================================================== the 12.7 design =======
fprintf('\n  ===== Worked design: 50 um water-in-oil droplets at 100 Hz =====\n');

w      = 50e-6;         % main channel width
hCh    = 50e-6;         % channel depth (for flow-rate bookkeeping only)
mu_c   = 10e-3;         % oil viscosity  [Pa.s]
sigma  = 5e-3;          % interfacial tension with surfactant [N/m]
ratio  = 0.5;           % Qd/Qc
fTarget= 100;           % Hz

alpha  = 1.0;
Lplug  = w*(1 + alpha*ratio);
Vdrop  = Lplug*w*hCh;                       % rectangular plug approximation

Qd = fTarget*Vdrop;                         % m^3/s
Qc = Qd/ratio;
Uc = Qc/(w*hCh);
Ca = mu_c*Uc/sigma;

fprintf('\n  Step 1-4 (geometry and flow ratio):\n');
fprintf('    w = %g um, Qd/Qc = %.2f -> L/w = %.2f, plug length %.1f um\n', ...
        w*1e6, ratio, 1+alpha*ratio, Lplug*1e6);
fprintf('    droplet volume ~ %.0f pL\n', Vdrop*1e15);
fprintf('\n  Step 5 (flow rates for %g Hz):\n', fTarget);
fprintf('    Qd = %.2f uL/min, Qc = %.2f uL/min\n', Qd*6e10, Qc*6e10);
fprintf('\n  Step 6 (THE CHECK):\n');
fprintf('    Uc = Qc/(w*h) = %.1f mm/s  ->  Ca = mu_c*Uc/sigma = %.4f\n', Uc*1e3, Ca);
if Ca < 1e-2
    fprintf('    Ca < 1e-2: SQUEEZING. Design is self-consistent.\n');
else
    fprintf('    *** Ca > 1e-2: this is DRIPPING, not squeezing. ***\n');
    fprintf('    The design is SELF-INCONSISTENT: we assumed the Garstecki\n');
    fprintf('    squeezing law to size the droplet, but the resulting flow\n');
    fprintf('    rate puts us outside its validity. Monodispersity will be\n');
    fprintf('    worse than specified.\n');
end

%% ---------------------------------------------- find a consistent design --
fprintf('\n  Step 7 (iterate): what frequency IS achievable in squeezing?\n');
CaMax  = 1e-2;
UcMax  = CaMax*sigma/mu_c;
QcMax  = UcMax*w*hCh;
QdMax  = ratio*QcMax;
fMax   = QdMax/Vdrop;
fprintf('    Ca <= %.0e  ->  Uc <= %.2f mm/s  ->  Qc <= %.2f uL/min\n', ...
        CaMax, UcMax*1e3, QcMax*6e10);
fprintf('    max frequency in squeezing = %.1f Hz  (target was %g Hz)\n', fMax, fTarget);
fprintf('\n    Options to recover the %g Hz target:\n', fTarget);
fprintf('      (a) less viscous oil: need mu_c <= %.1f mPa.s\n', ...
        mu_c*fMax/fTarget*1e3);
fprintf('      (b) more surfactant is the WRONG way -- lowering sigma RAISES Ca\n');
fprintf('      (c) parallelize: %d identical junctions on one chip\n', ...
        ceil(fTarget/fMax));
fprintf('      (d) accept dripping and a wider size distribution\n');
fprintf('    (c) is what real high-throughput devices do.\n');

%% ================================================== the regime table =====
fprintf('\n\n  ===== Regime map: Ca sweep at fixed geometry =====\n');
fprintf('  (mu_c = %g mPa.s, sigma = %g mN/m, w = %g um)\n\n', ...
        mu_c*1e3, sigma*1e3, w*1e6);
fprintf('     Uc [mm/s]        Ca        regime        predicted size\n');
fprintf('   -------------------------------------------------------------\n');
for U = [0.1 0.5 1 5 10 50 100]*1e-3
    C = mu_c*U/sigma;
    if     C < 1e-2, reg = 'squeezing'; sz = sprintf('L/w = %.2f (Garstecki)', 1+alpha*ratio);
    elseif C < 1e-1, reg = 'dripping '; sz = sprintf('d/w ~ %.2f (Ca^-1/3)', 0.5*C^(-1/3));
    else,            reg = 'jetting  '; sz = 'polydisperse';
    end
    fprintf('   %10.2f  %9.4f   %s   %s\n', U*1e3, C, reg, sz);
end

%% ================================== what would the simulation cost? ======
fprintf('\n\n  ===== Simulation cost (run this BEFORE you launch) =====\n');
rho1 = 998; rho2 = 900; nu = mu_c/rho2;
fprintf('\n   cells/w    dx [um]     dt_cap [s]    dt_visc [s]   steps for 3 drops\n');
fprintf('  ---------------------------------------------------------------------\n');
for nc = [10 20 40 80]
    dx = w/nc;
    dt_cap  = sqrt((rho1+rho2)*dx^3/(4*pi*sigma));
    dt_visc = 0.5/nu/(2/dx^2);
    dt = 0.4*min(dt_cap, dt_visc);
    tSim = 3/max(fMax,eps);
    fprintf('   %6d   %8.3f   %12.2e   %12.2e   %14.2e\n', ...
            nc, dx*1e6, dt_cap, dt_visc, tSim/dt);
end
fprintf(['\n   Tutorial 12.4 asks for >= 20 cells across the NECK, and the neck\n' ...
         '   is much thinner than w -- so the practical requirement is nearer\n' ...
         '   the bottom of this table than the top. Read off the step count and\n' ...
         '   decide whether you can afford it before you start.\n\n']);

if ~QUIET
    U = logspace(-5,-0.5,300);
    Ca_v = mu_c*U/sigma;
    figure('Color','w','Position',[80 80 640 420]);
    semilogx(Ca_v, ones(size(Ca_v)),'w'); hold on
    patch([1e-6 1e-2 1e-2 1e-6],[0 0 1 1],[0.3 0.7 0.4],'FaceAlpha',0.25,'EdgeColor','none');
    patch([1e-2 1e-1 1e-1 1e-2],[0 0 1 1],[0.9 0.8 0.3],'FaceAlpha',0.25,'EdgeColor','none');
    patch([1e-1 1e1  1e1  1e-1],[0 0 1 1],[0.9 0.4 0.3],'FaceAlpha',0.25,'EdgeColor','none');
    set(gca,'XScale','log'); xlim([1e-5 1]); ylim([0 1]);
    xline(Ca,'k--','LineWidth',2);
    text(3e-4,0.85,'SQUEEZING'); text(2.2e-2,0.85,'DRIP'); text(2e-1,0.85,'JETTING');
    text(Ca,0.15,sprintf('  your design Ca = %.3f',Ca));
    xlabel('Ca = \mu_c U_c / \sigma'); set(gca,'YTick',[]);
    title('T-junction regime map');
end
