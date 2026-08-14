%MIXING_STUDY  How far must the channel be for two streams to mix?
%
%   Tutorial 9.  Run:  mixing_study
%
%   Transports a passive scalar on a FROZEN flow field.  At Re << 1 the flow is
%   steady, so solving Navier-Stokes once and freezing it is exact, not an
%   approximation -- and it buys back the factor of Sc ~ 2000 by which the
%   viscous timestep limit is tighter than the scalar one (Tutorial 7.5).
%
%   Compares a straight channel with a serpentine of the same length, and --
%   critically -- MEASURES the numerical diffusivity of the grid actually used.
%   A mixing length is meaningless until you know the mesh is not producing it.
%
%   See also SCALAR_STEP, SERPENTINE, MIXING_INDEX, LIMITER.

if ~exist('QUIET','var'), QUIET = false; end
if ~exist('SCHEME','var'), SCHEME = 'vanleer'; end     % try 'upwind' for Ex 9.2

D     = 5e-10;
Umean = 1e-3;

%% ================================================== 1. straight channel ===
h  = 100e-6;  Lx = 4e-3;
nx = 400;  ny = 40;
dx = Lx/nx;  dy = h/ny;
xc = (dx/2:dx:Lx-dx/2)';  yc = (dy/2:dy:h-dy/2)';

uProf = 6*Umean*(yc/h).*(1-yc/h);
uS = repmat(uProf.', nx+1, 1);
vS = zeros(nx, ny+1);

Pe = Umean*h/D;
fprintf('\n  ============ MIXING STUDY ============\n');
fprintf('  Pe = %.0f, D = %.1e m^2/s, scheme = ''%s''\n', Pe, D, SCHEME);

%% ------------------------- numerical diffusion of THIS grid, measured -----
% The check that makes everything else believable.  Advect a Gaussian across
% the channel with D = 0 and see how much it spreads anyway.
sig0 = 6*dy;
cT = exp(-((yc-h/2).^2)/(2*sig0^2));
cTest = repmat(cT.', nx, 1).';  cTest = cTest.';        % nx x ny
uT = zeros(nx+1,ny);  vT = Umean*ones(nx,ny+1);         % advect in +y
dtT = 0.4*dy/Umean;
nT  = 40;
cc = cTest;
for n = 1:nT
    cc = scalar_step(cc, uT, vT, 0, dtT, dx, dy, SCHEME);
end
w  = max(cc(1,:),0).'; w = w/sum(w);
yb = mod(yc - sum(w.*yc) + h/2, h) - h/2;
sig1 = sqrt(sum(w.*yb.^2));
Dnum = max((sig1^2 - sig0^2)/(2*nT*dtT), 0);

fprintf('\n  Measured numerical diffusivity of this grid/scheme:\n');
fprintf('    dx = %.2f um, Pe_cell = %.2f\n', dx*1e6, Umean*dx/D);
fprintf('    D_num = %.3e m^2/s   D_num/D_phys = %.3f\n', Dnum, Dnum/D);
if Dnum/D > 0.1
    fprintf('    *** WARNING: above 0.1 -- the mixing result below is a\n');
    fprintf('    *** property of the MESH, not of the device. Refine, or\n');
    fprintf('    *** switch to a TVD limiter.\n');
else
    fprintf('    OK (below 0.1): physical diffusion dominates.\n');
end
fprintf('    For comparison, first-order upwind would give ~%.3e (%.1fx D).\n', ...
        Umean*dx/2, (Umean*dx/2)/D);

%% ------------------------------------------------- transport the scalar ---
c0S = double(repmat(yc.' > h/2, nx, 1));       % top half seeded
[MS, xS, ifaceS] = run_transport(c0S, uS, vS, D, dx, dy, xc, SCHEME, [], Umean);

%% ====================================================== 2. serpentine ====
flowFile = fullfile(fileparts(mfilename('fullpath')),'..','ch08','serpentine_flow.mat');
haveSerp = isfile(flowFile);
if haveSerp
    Sf = load(flowFile);
    Gs = Sf.G;  Ss = Sf.S;  solid = Sf.solid;
    ycs = Gs.yp;
    yCen = Sf.yCen;  halfw = Sf.halfw;

    % Seed the upper half of the channel cross-section at the inlet.
    c0 = zeros(Gs.nx, Gs.ny);
    for i = 1:Gs.nx
        yl = yCen(Gs.xp(i));
        c0(i,:) = double(ycs.' > yl & ~solid(i,:));
    end
    c0(solid) = 0;

    [MSerp, xSerp, ifaceSerp] = run_transport(c0, Ss.u, Ss.v, D, Gs.dx, Gs.dy, ...
                                              Gs.xp, SCHEME, solid, Sf.Umean);
else
    fprintf('\n  (serpentine_flow.mat not found -- run `serpentine` first to\n');
    fprintf('   include the serpentine comparison.)\n');
end

%% ---------------------------------------------------------- report -------
fprintf('\n  Mixing index M down the channel (M = 1 unmixed, 0 perfect):\n');
fprintf('    x [mm]  :'); fprintf(' %7.2f', xS(1:8:end)*1e3); fprintf('\n');
fprintf('    straight:'); fprintf(' %7.3f', MS(1:8:end));     fprintf('\n');
if haveSerp
    fprintf('    serpent.:'); fprintf(' %7.3f', MSerp(1:8:end)); fprintf('\n');
end

L95s = mix_length(xS, MS, 0.05);
fprintf('\n  Distance to M < 0.05:\n');
fprintf('    straight   : %s\n', fmt_len(L95s));
if haveSerp
    L95p = mix_length(xSerp, MSerp, 0.05);
    fprintf('    serpentine : %s\n', fmt_len(L95p));
end
fprintf('    0.3*Pe*w rule (Tutorial 4): %.2f mm\n', 0.3*Pe*h*1e3);

fprintf('\n  Interfacial length (relative to inlet), a leading indicator that\n');
fprintf('  the mixing index alone cannot see:\n');
fprintf('    straight   final/initial = %.3f\n', ifaceS(end)/ifaceS(1));
if haveSerp
    fprintf('    serpentine final/initial = %.3f\n', ifaceSerp(end)/ifaceSerp(1));
    fprintf(['\n  The serpentine improvement in 2-D is MODEST. That is the honest\n' ...
             '  result, not a bug: steady 2-D flow cannot be chaotic (streamlines\n' ...
             '  are level sets of the streamfunction, and particles cannot cross\n' ...
             '  them). Real mixers beat this using 3-D transverse flow --\n' ...
             '  herringbone grooves, Dean vortices -- which this model cannot\n' ...
             '  represent. See Exercise 9.6.\n']);
end
fprintf('\n');

if ~QUIET
    figure('Color','w','Position',[80 80 900 400]);
    plot(xS*1e3, MS,'LineWidth',1.8); hold on
    if haveSerp, plot(xSerp*1e3, MSerp,'LineWidth',1.8); end
    yline(0.05,'k--','M = 0.05');
    grid on; xlabel('downstream distance [mm]'); ylabel('mixing index M');
    if haveSerp, legend({'straight','serpentine'}); else, legend({'straight'}); end
    title(sprintf('Mixing, Pe = %.0f, scheme = %s', Pe, SCHEME));
end

%% ====================================================== local functions ===
function [M, xOut, iface] = run_transport(c0, u, v, D, dx, dy, xc, lim, solid, Umean)
%RUN_TRANSPORT  March the scalar to steady state, then read M(x).
nSteps = round(3.0*(xc(end)-xc(1))/Umean / (0.4*dx/max(abs(u(:)))));
dt = 0.4*min(dx/max(abs(u(:))), 1/(2*D)/(1/dx^2+1/dy^2));
c  = c0;
for n = 1:nSteps
    c = scalar_step(c, u, v, D, dt, dx, dy, lim, solid);
    c(1,:) = c0(1,:);                 % hold the inlet cross-section fixed
    c(end,:) = c(end-1,:);            % zero-gradient outlet
    if ~isempty(solid), c(solid) = 0; end
end

nx = size(c,1);
M  = zeros(nx,1);  iface = zeros(nx,1);
ref = c0(1,:);
for i = 1:nx
    if ~isempty(solid), open = ~solid(i,:); else, open = true(1,size(c,2)); end
    ci = c(i,open);
    if numel(ci) > 1 && std(ref(:)) > 0
        M(i) = std(ci)/std(ref(:));
    else
        M(i) = NaN;
    end
    iface(i) = sum(abs(diff(ci)));    % proxy for interfacial length
end
xOut = xc;
end

function L = mix_length(x, M, tol)
k = find(M < tol, 1);
if isempty(k), L = NaN; else, L = x(k); end
end

function s = fmt_len(L)
if isnan(L)
    s = 'not reached within the simulated length';
else
    s = sprintf('%.2f mm', L*1e3);
end
end
