%CAVITY  Lid-driven cavity: the standard first run of any NS solver.
%
%   Tutorial 7.  Run:  cavity
%
%   Runs the same geometry at an INERTIAL Reynolds number (100) and at a
%   MICROFLUIDIC one (0.1), so you can see what Re actually does to the flow.
%
%   See also NS_SOLVER, NS_DEFAULTS, CAVITY_VALIDATE.

if ~exist('QUIET','var'), QUIET = false; end

cases = struct('name', {'Re = 100 (inertial)', 'Re = 0.1 (microfluidic)'}, ...
               'Re',   {100, 0.1});
out = cell(1,2);

for k = 1:2
    opt = ns_defaults();
    opt.nx = 64;  opt.ny = 64;
    opt.Lx = 1;   opt.Ly = 1;
    opt.rho = 1;  opt.nu = 1/cases(k).Re;
    opt.bc.top = struct('type','moving','u',1);
    opt.tEnd   = 60;
    opt.steady = true;  opt.steadyTol = 2e-5;
    opt.report = 500;

    fprintf('\n  ========== %s ==========\n', cases(k).name);
    out{k} = ns_solver(opt);

    S = out{k};  G = S.G;
    % Locate the primary vortex from the streamfunction extremum.
    psi = streamfun(S);
    [~,idx] = max(abs(psi(:)));
    [i0,j0] = ind2sub(size(psi), idx);
    fprintf('    primary vortex centre: (x,y) = (%.3f, %.3f)\n', G.xp(i0), G.yp(j0));
    fprintf('    |psi|_max = %.4e\n', abs(psi(idx)));
end

fprintf(['\n  At Re = 100 the vortex centre sits DOWNSTREAM of and below the\n' ...
         '  geometric centre -- inertia carries the fluid past the corner.\n' ...
         '  At Re = 0.1 it moves back toward the centre and the flow becomes\n' ...
         '  nearly fore-aft symmetric. That symmetry IS the reversibility of\n' ...
         '  Stokes flow from Tutorial 1, visible directly in the streamlines.\n\n']);

if ~QUIET
    figure('Color','w','Position',[80 80 940 440]);
    for k = 1:2
        S = out{k};  G = S.G;
        subplot(1,2,k);
        [Xc,Yc] = ndgrid(G.xp, G.yp);
        streamslice(Xc.', Yc.', S.uc.', S.vc.', 3);
        axis equal tight; xlabel('x'); ylabel('y');
        title(cases(k).name);
    end
end
