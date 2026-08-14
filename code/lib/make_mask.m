function m = make_mask(G, varargin)
%MAKE_MASK  Build a solid-cell mask from geometric primitives.
%
%   m = make_mask(G, shape1, shape2, ...)   returns an nx-by-ny logical array,
%   true in cells that are SOLID.  Shapes are unioned together.
%
%   Each shape is a cell array:
%       {'rect',   x0, x1, y0, y1}     axis-aligned rectangle
%       {'circle', xc, yc, r}          disc
%       {'below',  @(x) yfun(x)}       everything with y < yfun(x)
%       {'above',  @(x) yfun(x)}       everything with y > yfun(x)
%       {'fun',    @(X,Y) logicalExpr} arbitrary predicate on ndgrid arrays
%
%   Coordinates are physical (metres), evaluated at CELL CENTRES.
%
%   Example -- a channel with a cylindrical post:
%       G = mac_grid(200,50,2e-3,5e-4);
%       m = make_mask(G, {'circle', 5e-4, 2.5e-4, 8e-5});
%
%   Example -- a serpentine, as walls above and below a wavy centreline:
%       amp = 1e-4; lam = 5e-4; halfw = 6e-5;
%       yc  = @(x) 2.5e-4 + amp*sin(2*pi*x/lam);
%       m   = make_mask(G, {'below', @(x) yc(x)-halfw}, ...
%                          {'above', @(x) yc(x)+halfw});
%
%   RESOLUTION WARNING.  The mask is staircased: a circle becomes a pixelated
%   blob.  Use at least 20 cells across the smallest feature, and remember that
%   forces computed on a staircased body are first-order accurate at best (see
%   Tutorial 8.3).  make_mask warns if a 'circle' is thinner than 10 cells.
%
%   See also MAC_GRID, NS_SOLVER, OBSTACLE, SERPENTINE.

X = G.Xp;  Y = G.Yp;                 % cell centres, ndgrid order (x,y)
m = false(G.nx, G.ny);

for k = 1:numel(varargin)
    s = varargin{k};
    switch lower(s{1})
        case 'rect'
            m = m | (X >= s{2} & X <= s{3} & Y >= s{4} & Y <= s{5});

        case 'circle'
            r = s{4};
            if 2*r < 10*min(G.dx,G.dy)
                warning('make_mask:coarse', ...
                    ['Circle of radius %.3g m spans only %.1f cells. ' ...
                     'Use >= 20 across the smallest feature.'], ...
                    r, 2*r/min(G.dx,G.dy));
            end
            m = m | ((X-s{2}).^2 + (Y-s{3}).^2 <= r^2);

        case 'below'
            m = m | (Y < s{2}(X));

        case 'above'
            m = m | (Y > s{2}(X));

        case 'fun'
            m = m | logical(s{2}(X,Y));

        otherwise
            error('make_mask:shape','Unknown shape ''%s''.', s{1});
    end
end
end
