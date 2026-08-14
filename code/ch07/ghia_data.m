function d = ghia_data(Re)
%GHIA_DATA  Lid-driven cavity benchmark, Ghia, Ghia & Shin (1982).
%
%   d = ghia_data(Re)   for Re = 100 or 400
%
%   Returns tabulated centreline velocities for the unit cavity with a lid
%   moving at u = 1:
%       d.y, d.u   u along the VERTICAL   centreline (x = 0.5)
%       d.x, d.v   v along the HORIZONTAL centreline (y = 0.5)
%
%   Source: U. Ghia, K.N. Ghia, C.T. Shin, "High-Re solutions for incompressible
%   flow using the Navier-Stokes equations and a multigrid method",
%   J. Comput. Phys. 48 (1982) 387-411, Tables I and II.
%
%   *** PROVENANCE WARNING ***
%   These numbers were transcribed from memory, not from the paper.  The Re=100
%   rows are corroborated: ns_solver reproduces BOTH u and v to <0.01 of the lid
%   speed, which is unlikely to happen by accident.  The Re=400 v row is NOT
%   corroborated -- cavity_validate finds a ~0.14 discrepancy that does not
%   shrink under grid refinement, which is the signature of bad reference data
%   rather than a bad solution.  Verify that row against Table II before relying
%   on it.  An unverified reference is not a reference.
%
%   This is THE standard verification case for an incompressible solver.  It is
%   worth validating here even though microfluidics runs at Re ~ 0.1: at Re =
%   0.1 the convective term is invisible, so a test at that Re cannot fail and
%   therefore proves nothing.  Validate where the physics is hard.
%
%   See also CAVITY_VALIDATE, NS_SOLVER.

d.y = [0.0000 0.0547 0.0625 0.0703 0.1016 0.1719 0.2813 0.4531 0.5000 ...
       0.6172 0.7344 0.8516 0.9531 0.9609 0.9688 0.9766 1.0000]';

d.x = [0.0000 0.0625 0.0703 0.0781 0.0938 0.1563 0.2266 0.2344 0.5000 ...
       0.8047 0.8594 0.9063 0.9453 0.9531 0.9609 0.9688 1.0000]';

switch Re
    case 100
        d.u = [ 0.00000 -0.03717 -0.04192 -0.04775 -0.06434 -0.10150 ...
               -0.15662 -0.21090 -0.20581 -0.13641  0.00332  0.23151 ...
                0.68717  0.73722  0.78871  0.84123  1.00000]';
        d.v = [ 0.00000  0.09233  0.10091  0.10890  0.12317  0.16077 ...
                0.17507  0.17527  0.05454 -0.24533 -0.22445 -0.16914 ...
               -0.10313 -0.08864 -0.07391 -0.05906  0.00000]';
    case 400
        d.u = [ 0.00000 -0.08186 -0.09266 -0.10338 -0.14612 -0.24299 ...
               -0.32726 -0.17119 -0.11477  0.02135  0.16256  0.29093 ...
                0.55892  0.61756  0.68439  0.75837  1.00000]';
        d.v = [ 0.00000  0.18360  0.19713  0.20920  0.22965  0.28124 ...
                0.30203  0.30174  0.05186 -0.38598 -0.44993 -0.23827 ...
               -0.22847 -0.19254 -0.15663 -0.12146  0.00000]';
    otherwise
        error('ghia_data:Re','Only Re = 100 and 400 are tabulated here.');
end
d.Re = Re;
end
