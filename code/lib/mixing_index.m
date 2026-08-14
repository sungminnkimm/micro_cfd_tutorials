function M = mixing_index(c, cRef)
%MIXING_INDEX  Normalized concentration standard deviation.
%
%   M = mixing_index(c, cRef)
%
%   M = std(c) / std(cRef), the standard measure of how well two streams have
%   mixed.  With cRef the initial (fully segregated) field:
%
%       M = 1   completely unmixed
%       M = 0   perfectly uniform
%
%   Convention note: some papers report the MIXING EFFICIENCY 1 - M instead, so
%   a "95% mixed" claim can mean M = 0.05 or M = 0.95 depending on the author.
%   State which one you mean.  This course always reports M as defined above.
%
%   Caveat: M is blind to structure.  A field striped into many thin lamellae
%   of pure A and pure B has the same std as one thick stripe of each, but is
%   about to mix far faster.  When you stretch and fold (Tutorial 9), watch the
%   interfacial length as well -- M alone will understate what the mixer did
%   until diffusion finally cashes in the extra interface.
%
%   For a weighted grid (non-uniform cells), pass volume-weighted fields or use
%   a weighted std; on the uniform grids in this course, plain std is correct.
%
%   See also DIFFUSION1D, MIXING_STUDY.

if nargin < 2 || isempty(cRef)
    error('mixing_index:ref','Provide a reference field (usually the initial one).');
end

s0 = std(cRef(:));
if s0 == 0
    error('mixing_index:flat','Reference field is uniform; std is 0.');
end
M = std(c(:)) / s0;
end
