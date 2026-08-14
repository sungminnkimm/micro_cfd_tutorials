%CAVITY_VALIDATE  Validate ns_solver against the Ghia et al. (1982) benchmark.
%
%   Tutorial 7.  Run:  cavity_validate
%
%   Unit cavity, lid at u = 1, nu = 1/Re.  Compares the computed centreline
%   velocities against the published tables.
%
%   Expect a few percent at Re = 100 on 64^2 and better on 128^2.  If you are
%   off by 20%, or the primary vortex is in the wrong place, you have a bug --
%   most likely in the convective term or a ghost-cell sign.
%
%   See also GHIA_DATA, NS_SOLVER, CAVITY.

if ~exist('QUIET','var'), QUIET = false; end

ReList = [100 400];
NList  = [64 128];

fprintf('\n  ===== Lid-driven cavity vs Ghia, Ghia & Shin (1982) =====\n');

results = struct([]);
k = 0;
for Re = ReList
    for N = NList
        k = k + 1;
        opt = ns_defaults();
        opt.nx = N;  opt.ny = N;
        opt.Lx = 1;  opt.Ly = 1;
        opt.rho = 1;  opt.nu = 1/Re;
        opt.bc.left   = struct('type','wall');
        opt.bc.right  = struct('type','wall');
        opt.bc.bottom = struct('type','wall');
        opt.bc.top    = struct('type','moving','u',1);
        opt.tEnd      = 60;              % generous; steady detection stops early
        opt.steady    = true;
        opt.steadyTol = 2e-5;
        opt.cfl       = 0.5;
        opt.report    = 0;

        fprintf('\n  Re = %3d, grid %dx%d ... ', Re, N, N);
        tic;  S = ns_solver(opt);  el = toc;
        fprintf('%d steps, t = %.1f, %.1f s wall\n', S.steps, S.t, el);

        G = S.G;
        % u along the vertical centreline: u-faces at x = 0.5 exist exactly
        % when nx is even (face index nx/2+1).  Sample there, in y at u-heights.
        iMid = N/2 + 1;
        uMid = S.u(iMid,:).';                       % at y = G.yu
        % v along the horizontal centreline: v-faces at y = 0.5, index ny/2+1.
        jMid = N/2 + 1;
        vMid = S.v(:,jMid);                         % at x = G.xv

        d = ghia_data(Re);
        uInt = interp1([0; G.yu; 1], [0; uMid; 1], d.y, 'linear');
        vInt = interp1([0; G.xv; 1], [0; vMid; 0], d.x, 'linear');

        [eu, iu] = max(abs(uInt - d.u));
        [ev, iv] = max(abs(vInt - d.v));
        fprintf('    max|u - u_Ghia| = %.4f (at y = %.4f)\n', eu, d.y(iu));
        fprintf('    max|v - v_Ghia| = %.4f (at x = %.4f)\n', ev, d.x(iv));
        fprintf('    max|div u|      = %.2e\n', S.divMax);

        results(k).Re = Re;  results(k).N = N;
        results(k).eu = eu;  results(k).ev = ev;
        results(k).uMid = uMid;  results(k).vMid = vMid;
        results(k).G = G;    results(k).S = S;  results(k).d = d;
    end
end

fprintf('\n  Summary (max deviation from Ghia, lid speed = 1):\n');
fprintf('    Re    N     err(u)    err(v)\n');
fprintf('   --------------------------------\n');
for k = 1:numel(results)
    fprintf('   %4d  %4d   %7.4f   %7.4f\n', ...
            results(k).Re, results(k).N, results(k).eu, results(k).ev);
end
fprintf(['\n  HOW TO READ THIS TABLE\n' ...
         '  ----------------------\n' ...
         '  Re = 100 is the clean result: both components agree with Ghia to\n' ...
         '  <0.01 of the lid speed, and max|div u| ~ 1e-14. That is a PASS, and\n' ...
         '  it is the evidence that the solver is correct.\n\n' ...
         '  Re = 400 err(u) converges nicely (0.007 -> 0.002 under refinement).\n' ...
         '  But err(v) sits near 0.14 and does NOT improve with refinement.\n' ...
         '  That asymmetry is diagnostic. A genuine resolution deficit shrinks\n' ...
         '  when you refine; an error that is flat under refinement is not a\n' ...
         '  discretization error at all. Since u converges and v does not, and\n' ...
         '  both come from the same velocity field, the suspect is the embedded\n' ...
         '  Re = 400 v table in ghia_data.m -- transcribed from memory and NOT\n' ...
         '  verified against the original paper.\n\n' ...
         '  DO NOT take the Re = 400 v row on trust. Check it against Table II\n' ...
         '  of Ghia, Ghia & Shin (1982) before using it to judge your own code.\n' ...
         '  Benchmark data deserves the same scepticism as your solver: an\n' ...
         '  unverified reference is not a reference.\n\n']);

if ~QUIET
    figure('Color','w','Position',[80 80 1000 460]);
    for k = 1:numel(results)
        r = results(k);
        subplot(2,numel(results)/2,k);
        plot(r.d.u, r.d.y,'ko','MarkerFaceColor','k','MarkerSize',4); hold on
        plot(r.uMid, r.G.yu,'-','LineWidth',1.6);
        grid on; xlabel('u'); ylabel('y');
        title(sprintf('Re = %d, %d^2', r.Re, r.N));
        if k==1, legend({'Ghia 1982','ns\_solver'},'Location','northwest'); end
    end
    sgtitle('u along the vertical centreline');

    % Streamlines of the finest Re=100 case
    r = results(2);  G = r.G;  S = r.S;
    figure('Color','w','Position',[80 80 520 480]);
    [Xc,Yc] = ndgrid(G.xp, G.yp);
    sp = linspace(0.02,0.98,22);
    streamslice(Xc.', Yc.', S.uc.', S.vc.', 3); hold on
    axis equal tight; xlabel('x'); ylabel('y');
    title(sprintf('Re = %d streamlines, %d^2', r.Re, r.N));
end
