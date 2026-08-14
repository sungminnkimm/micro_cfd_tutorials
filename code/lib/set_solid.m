function opt = set_solid(opt, solid)
%SET_SOLID  Attach a solid mask and its face masks to an options struct.
%
%   opt = set_solid(opt, solid)
%
%   Expands the nx-by-ny cell mask onto the staggered faces: a face is blocked
%   if EITHER neighbouring cell is solid.  Used by NS2PHASE (NS_SOLVER does the
%   same expansion internally).
%
%   See also MAKE_MASK, NS2PHASE, NS_SOLVER.

nx = opt.nx;  ny = opt.ny;
if isempty(solid)
    opt.solid = [];  opt.uMask = [];  opt.vMask = [];
    return
end
if ~isequal(size(solid),[nx ny])
    error('set_solid:size','solid must be %dx%d.', nx, ny);
end
solid = logical(solid);

uMask = false(nx+1, ny);
uMask(1:end-1,:) = uMask(1:end-1,:) | solid;
uMask(2:end,:)   = uMask(2:end,:)   | solid;

vMask = false(nx, ny+1);
vMask(:,1:end-1) = vMask(:,1:end-1) | solid;
vMask(:,2:end)   = vMask(:,2:end)   | solid;

opt.solid = solid;  opt.uMask = uMask;  opt.vMask = vMask;
end
