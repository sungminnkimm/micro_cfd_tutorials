%RUN_ALL_TESTS  Execute every course script and report pass/fail.
%
%   cd D:\Projects\CFD\code
%   run_all_tests
%
%   Runs each tutorial's scripts with plotting suppressed and checks the key
%   numeric claims.  Use it as a regression test after you modify anything in
%   lib/.
%
%   Note the deliberately long-running entries at the end; pass 'fast' to skip
%   them:   run_all_tests fast

if ~exist('MODE','var'), MODE = 'full'; end

setup_path;
QUIET = true;                                    %#ok<NASGU> used by the scripts
tAll = tic;

tests = {
  'ch01 scaling'         'scaling_demo'
  'ch01 sweep'           'scale_sweep'
  'ch02 convergence'     'convergence_study'
  'ch02 stencils'        'fd_stencils'
  'ch03 Poisson MMS'     'poisson_mms'
  'ch03 conduction'      'poisson_demo'
  'ch04 stability'       'stability_demo'
  'ch04 diffusion 1-D'   'diffusion1d'
  'ch04 diffusion 2-D'   'diffusion2d'
  'ch05 advection'       'advect1d'
  'ch05 cell Peclet'     'cell_peclet'
  'ch05 num diffusion'   'num_diffusion'
  'ch06 checkerboard'    'checkerboard_demo'
  'ch06 Poiseuille'      'poiseuille'
  'ch07 cavity'          'cavity'
  'ch08 channel'         'channel'
  'ch09 Taylor-Aris'     'taylor_dispersion'
  'ch10 Young-Laplace'   'young_laplace'
  'ch10 capillary rise'  'capillary_rise'
  'ch10 meniscus'        'meniscus'
  'ch11 static droplet'  'static_droplet'
  'ch12 regime map'      'regime_map'
};

% Long-running; skipped by  run_all_tests fast
slow = {
  'ch07 Ghia validation' 'cavity_validate'
  'ch08 obstacle'        'obstacle'
  'ch08 serpentine'      'serpentine'
  'ch09 mixing study'    'mixing_study'
  'ch11 rising bubble'   'rising_bubble'
  'ch12 T-junction'      't_junction'
  'ch12 flow focusing'   'flow_focusing'
  'ch13 MMS driver'      'mms_driver'
};

if ~strcmpi(MODE,'fast')
    tests = [tests; slow];
end

nPass = 0;  nFail = 0;  failed = {};
fprintf('\n================ RUNNING %d TESTS ================\n', size(tests,1));

for k = 1:size(tests,1)
    name = tests{k,1};  script = tests{k,2};
    fprintf('\n--- [%2d/%2d] %-22s (%s)\n', k, size(tests,1), name, script);
    t0 = tic;
    try
        evalin('base', sprintf('QUIET=true; %s;', script));
        fprintf('    PASS  (%.1f s)\n', toc(t0));
        nPass = nPass + 1;
    catch ME
        fprintf(2, '    FAIL  %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf(2, '          at %s line %d\n', ME.stack(1).name, ME.stack(1).line);
        end
        nFail = nFail + 1;  failed{end+1} = name; %#ok<SAGROW>
    end
end

fprintf('\n================ SUMMARY ================\n');
fprintf('  %d passed, %d failed, total %.1f s\n', nPass, nFail, toc(tAll));
if nFail > 0
    fprintf('  failures:\n');
    fprintf('    %s\n', failed{:});
else
    fprintf('  ALL SCRIPTS RAN CLEAN.\n');
end
fprintf('=========================================\n\n');
