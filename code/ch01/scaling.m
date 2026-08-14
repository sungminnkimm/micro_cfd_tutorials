function s = scaling(varargin)
%SCALING  Dimensionless groups and derived scales for a microscale flow.
%
%   s = scaling('L',100e-6, 'U',1e-3)                 % water defaults
%   s = scaling('L',50e-6,  'U',5e-3, 'D',1e-10)
%   s = scaling(...,'rho',1000,'mu',1e-3,'sigma',0.072,'g',9.81)
%
%   Returns a struct with fields Re, Pe, Ca, Bo, We, Sc, plus:
%     t_visc   momentum diffusion time across L       (rho*L^2/mu)
%     t_diff   mass diffusion time across L           (L^2/D)
%     t_adv    advection time across L                (L/U)
%     t_cap    capillary wave time on scale L         (sqrt(rho*L^3/sigma))
%     L_mix    channel length for diffusive mixing    (Pe*L)
%     l_cap    capillary length                       (sqrt(sigma/(rho*g)))
%     d_ent    hydrodynamic entrance length           (0.06*Re*L)
%     dp_pois  pressure gradient for this U in a slot (12*mu*U/L^2)
%
%   See also SCALING_DEMO, SCALE_SWEEP.

p = inputParser;
p.addParameter('L',     100e-6);   % m      characteristic length
p.addParameter('U',     1e-3);     % m/s    characteristic velocity
p.addParameter('rho',   998);      % kg/m^3
p.addParameter('mu',    1.0e-3);   % Pa.s
p.addParameter('sigma', 0.072);    % N/m
p.addParameter('D',     5e-10);    % m^2/s  small molecule in water
p.addParameter('g',     9.81);     % m/s^2
p.parse(varargin{:});
q = p.Results;

L = q.L;  U = q.U;  rho = q.rho;  mu = q.mu;  sig = q.sigma;  D = q.D;  g = q.g;
nu = mu/rho;

s    = q;
s.nu = nu;

% --- dimensionless groups ------------------------------------------------
s.Re = rho*U*L/mu;          % inertia   / viscous
s.Pe = U*L/D;               % advection / mass diffusion
s.Ca = mu*U/sig;            % viscous   / surface tension
s.Bo = rho*g*L^2/sig;       % gravity   / surface tension
s.We = rho*U^2*L/sig;       % inertia   / surface tension  (= Re*Ca)
s.Sc = nu/D;                % momentum diffusivity / mass diffusivity

% --- time scales ---------------------------------------------------------
s.t_visc = L^2/nu;              % momentum diffuses across L
s.t_diff = L^2/D;               % molecules diffuse across L
s.t_adv  = L/U;                 % fluid crosses L
s.t_cap  = sqrt(rho*L^3/sig);   % capillary wave period on scale L

% --- length scales -------------------------------------------------------
s.L_mix = s.Pe*L;                    % channel length to mix by diffusion
s.l_cap = sqrt(sig/(rho*g));         % capillary length
s.d_ent = 0.06*s.Re*L;               % hydrodynamic entrance length (laminar)

% --- engineering ---------------------------------------------------------
s.dp_pois = 12*mu*U/L^2;             % dP/dx for mean U in a slot of height L
end
