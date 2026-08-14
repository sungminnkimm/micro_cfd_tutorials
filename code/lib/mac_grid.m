function G = mac_grid(nx, ny, Lx, Ly)
%MAC_GRID  Build a staggered (Marker-and-Cell) grid descriptor.
%
%   G = mac_grid(nx, ny, Lx, Ly)
%
%   Variable layout for nx-by-ny CELLS:
%
%       p : nx    x  ny       cell centres
%       u : (nx+1) x ny       vertical   faces, includes BOTH x-boundaries
%       v : nx    x (ny+1)    horizontal faces, includes BOTH y-boundaries
%
%   The "+1" is the entire source of staggered-grid bugs: a row of nx cells has
%   nx+1 faces (fenceposts).  Keep this picture in mind:
%
%              v(i,j+1)
%                 ^
%          +------|------+
%          |             |
%     u(i,j) -> p(i,j) -> u(i+1,j)
%          |             |
%          +------|------+
%                 ^
%              v(i,j)
%
%   WHICH BOUNDARIES NEED GHOST CELLS
%   ---------------------------------
%   u is stored ON the x-walls, so a no-slip condition at x=0 or x=Lx is just
%   u(1,:) = 0 or u(end,:) = 0 -- an actual unknown, set directly.  But u is
%   NOT stored on the y-walls (it lives at cell-centre height), so no-slip at
%   y=0 or y=Ly needs a ghost: u_ghost = 2*u_wall - u_adjacent.
%
%   For v it is exactly the other way round.  Half your boundaries are direct
%   and half need ghosts, and which is which SWAPS between u and v.  Getting
%   this backwards is the single most common bug in a new staggered solver.
%
%       variable   x=0, x=Lx      y=0, y=Ly
%       --------   ------------   ------------
%       u          direct         ghost needed
%       v          ghost needed   direct
%       p          ghost (BCs enter through the Poisson RHS, Tutorial 3)
%
%   FIELDS
%       nx ny Lx Ly dx dy
%       xp yp    cell-centre coordinates (nx, ny)
%       xu yu    u-face coordinates      (nx+1, ny)
%       xv yv    v-face coordinates      (nx, ny+1)
%       Xp Yp    ndgrid arrays at cell centres
%       Xu Yu    ndgrid arrays at u-faces
%       Xv Yv    ndgrid arrays at v-faces
%
%   See also POISEUILLE, NS_SOLVER, CHECKERBOARD_DEMO.

G.nx = nx;  G.ny = ny;  G.Lx = Lx;  G.Ly = Ly;
G.dx = Lx/nx;  G.dy = Ly/ny;

G.xp = (G.dx/2 : G.dx : Lx - G.dx/2)';    % nx
G.yp = (G.dy/2 : G.dy : Ly - G.dy/2)';    % ny
G.xu = (0      : G.dx : Lx)';             % nx+1
G.yu = G.yp;                              % ny    (u sits at cell-centre height)
G.xv = G.xp;                              % nx
G.yv = (0      : G.dy : Ly)';             % ny+1

[G.Xp, G.Yp] = ndgrid(G.xp, G.yp);
[G.Xu, G.Yu] = ndgrid(G.xu, G.yu);
[G.Xv, G.Yv] = ndgrid(G.xv, G.yv);

% Sizes, kept explicitly so callers can assert against them.
G.size_p = [nx   ny  ];
G.size_u = [nx+1 ny  ];
G.size_v = [nx   ny+1];
end
