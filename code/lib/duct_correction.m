function f = duct_correction(ar, nTerms)
%DUCT_CORRECTION  How much of a 2-D flow-rate prediction survives in 3-D.
%
%   f = duct_correction(w_over_h)
%
%   Returns  ubar_duct / ubar_2D  for pressure-driven laminar flow in a
%   rectangular duct of width w and height h, from the exact series solution
%
%       ubar = (G h^2)/(12 mu) * [ 1 - sum_{n odd} 192/(n^5 pi^5) * (h/w)
%                                        * tanh(n pi w / (2h)) ]
%
%   The bracket is f.  Multiply any 2-D flow rate or mean velocity by it to get
%   the rectangular-duct value; divide a required flow rate by it to size a
%   channel.  Equivalently, the pressure drop for a given flow rate is 1/f
%   times the 2-D prediction.
%
%   WHY YOU NEED THIS
%       ar = 1   -> f = 0.42   2-D overpredicts flow rate by 137%
%       ar = 5   -> f = 0.86               by 16%
%       ar = 10  -> f = 0.94               by 7%
%       ar -> Inf-> f = 1      2-D is exact
%
%   Every 2-D simulation in this course predicts flow MAGNITUDE too high for a
%   real channel, and this is the correction.  Flow STRUCTURE (recirculation,
%   streamline topology, where a droplet pinches) survives 2-D much better than
%   magnitude does.
%
%   ar may be an array.  nTerms defaults to 50 odd terms (converges fast; the
%   n^-5 weight means term 11 is already ~1e-5 of term 1).
%
%   See also POISEUILLE.

if nargin < 2, nTerms = 50; end

n = 1:2:(2*nTerms-1);                  % odd integers
f = zeros(size(ar));

for k = 1:numel(ar)
    a = ar(k);
    % tanh saturates to 1 for large n*a; that is fine and exact in double.
    s = sum( 192./(n.^5*pi^5) .* (1/a) .* tanh(n*pi*a/2) );
    f(k) = 1 - s;
end
end
