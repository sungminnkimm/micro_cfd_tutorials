function txt = ch01_verdict(val, lims, lowTxt, highTxt, midTxt)
%CH01_VERDICT  Small helper: turn a number and two thresholds into a comment.
%
%   ch01_verdict(0.1, [1 100], 'Stokes', 'inertial', 'transitional')
%       -> '<-- Stokes'

if     val < lims(1), txt = ['<-- ' lowTxt];
elseif val > lims(2), txt = ['<-- ' highTxt];
else,                 txt = ['<-- ' midTxt];
end
end
