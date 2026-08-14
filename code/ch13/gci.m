function out = gci(h, f, Fs)
%GCI  Grid Convergence Index: a numerical error bar you can publish.
%
%   out = gci(h, f)        h = [h1 h2 h3] grid spacings, FINEST FIRST
%                          f = [f1 f2 f3] the quantity of interest on each
%   out = gci(h, f, Fs)    safety factor (default 1.25 for 3 grids)
%
%   Richardson extrapolation with refinement ratio r = h2/h1:
%       p     = ln|(f3-f2)/(f2-f1)| / ln(r)          observed order
%       f_ext = f1 + (f1-f2)/(r^p - 1)               extrapolated exact value
%       GCI   = Fs*|f1-f2| / ((r^p-1)*|f1|) * 100%   fine-grid error bar
%
%   THE CHECK PEOPLE SKIP: the observed p must come out near your THEORETICAL
%   order.  If you get p = 0.4 or p = 6, the grids are not in the asymptotic
%   range, Richardson extrapolation is invalid, and the error bar it produces
%   is meaningless.  Refine further before believing it.  out.asymptotic flags
%   this.
%
%   Reporting "3.4 mm with a grid convergence index of 2.1%" is a defensible
%   claim.  "We used a fine mesh" is not.
%
%   Example:
%       out = gci([1 2 4]*1e-6, [3.41 3.52 3.88]);
%
%   See also MMS_DRIVER.

if nargin < 3, Fs = 1.25; end
h = h(:).';  f = f(:).';
if numel(h)~=3 || numel(f)~=3
    error('gci:size','Need exactly three grids, finest first.');
end
if h(1) >= h(2) || h(2) >= h(3)
    error('gci:order','h must be strictly increasing coarseness: h1 < h2 < h3.');
end

r21 = h(2)/h(1);
r32 = h(3)/h(2);

e21 = f(2)-f(1);
e32 = f(3)-f(2);

if abs(e21) < eps*max(abs(f))
    warning('gci:converged','f1 and f2 are identical to round-off; p undefined.');
    out.p = NaN;
else
    if abs(r21 - r32) > 1e-8
        % Non-constant refinement ratio needs the implicit Roache formula.
        pfun = @(p) p - abs(log(abs(e32/e21)) + log((r21^p-1)/(r32^p-1)))/log(r21);
        out.p = fzero(pfun, 2);
    else
        out.p = log(abs(e32/e21))/log(r21);
    end
end

p = out.p;
out.fExtrap = f(1) + e21/(r21^p - 1);
out.gci21   = Fs*abs(e21/f(1))/(r21^p - 1)*100;
out.gci32   = Fs*abs(e32/f(2))/(r32^p - 1)*100;
out.r       = r21;
out.f       = f;
out.h       = h;

% Asymptotic-range check: the two GCIs should be related by r^p.
out.asymRatio  = out.gci32/(r21^p * out.gci21);
out.asymptotic = abs(out.asymRatio - 1) < 0.1;

fprintf('\n  Grid Convergence Index\n');
fprintf('    h  = [%.4g %.4g %.4g]\n', h);
fprintf('    f  = [%.6g %.6g %.6g]\n', f);
fprintf('    observed order p     = %.3f\n', p);
fprintf('    extrapolated value   = %.6g\n', out.fExtrap);
fprintf('    GCI(fine)            = %.3f %%\n', out.gci21);
fprintf('    GCI(coarse)          = %.3f %%\n', out.gci32);
fprintf('    asymptotic ratio     = %.3f  (want ~1)\n', out.asymRatio);
if out.asymptotic
    fprintf('    -> in the asymptotic range: the error bar is meaningful.\n\n');
else
    fprintf(['    -> NOT in the asymptotic range. Richardson extrapolation is\n' ...
             '       invalid here and this error bar means nothing. Refine\n' ...
             '       further before quoting it.\n\n']);
end
end
