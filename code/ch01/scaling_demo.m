%SCALING_DEMO  Print the dimensionless picture for a typical microchannel.
%
%   Tutorial 1.  Run:  scaling_demo
%
%   See also SCALING, SCALE_SWEEP.

if ~exist('QUIET','var'), QUIET = false; end

s = scaling('L',100e-6,'U',1e-3);

fprintf('\n--- Microchannel: L = %g um, U = %g mm/s (water) ---\n', ...
        s.L*1e6, s.U*1e3);
fprintf('  Re = %10.3g   %s\n', s.Re, ch01_verdict(s.Re, [1 100], ...
        'Stokes flow: drop inertia', 'inertia matters', 'transitional'));
fprintf('  Pe = %10.3g   %s\n', s.Pe, ch01_verdict(s.Pe, [1 100], ...
        'diffusion-dominated', 'advection-dominated: mixing is HARD', 'comparable'));
fprintf('  Ca = %10.3g   %s\n', s.Ca, ch01_verdict(s.Ca, [1e-3 1e-1], ...
        'surface tension dominates interfaces', 'shear deforms interfaces', 'transitional'));
fprintf('  Bo = %10.3g   %s\n', s.Bo, ch01_verdict(s.Bo, [1e-2 1], ...
        'gravity negligible', 'gravity matters', 'comparable'));
fprintf('  We = %10.3g\n', s.We);
fprintf('  Sc = %10.3g   (momentum diffuses %.0fx faster than mass)\n', s.Sc, s.Sc);

fprintf('\n  Time scales:\n');
fprintf('    momentum diffusion across L : %10.3g s\n', s.t_visc);
fprintf('    advection across L          : %10.3g s\n', s.t_adv);
fprintf('    capillary wave              : %10.3g s\n', s.t_cap);
fprintf('    MASS diffusion across L     : %10.3g s   <-- the slow one\n', s.t_diff);

fprintf('\n  Design numbers:\n');
fprintf('    channel length to mix by diffusion : %8.3g mm\n', s.L_mix*1e3);
fprintf('    capillary length (water)           : %8.3g mm\n', s.l_cap*1e3);
fprintf('    entrance length                    : %8.3g um\n', s.d_ent*1e6);
fprintf('    dP/dx to drive this flow           : %8.3g bar/cm\n', s.dp_pois*1e-5*1e-2);
fprintf('\n');
