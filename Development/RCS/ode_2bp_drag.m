function dy = ode_2bp_drag(t, y, vars)
% ode_2bp_perturbed.m : Returns the ODE system for the Two-Body problem in 
% cartesian coordinates (state vector representation) in presence of drag
% perturbation.
%
% INPUT:
%   t           [1]         =   Time variable (required for ode MATLAB schemes to work) [s] 
%   y           [6x1]       =   State vector
%   y(1)                    =   Position along x direction w.r.t massive body           [km]
%   y(2)                    =   Position along y direction w.r.t massive body           [km]
%   y(3)                    =   Position along z direction w.r.t massive body           [km]
%   y(4)                    =   Velocity along x direction w.r.t massive body           [km/s]
%   y(5)                    =   Velocity along y direction w.r.t massive body           [km/s]
%   y(6)                    =   Velocity along z direction w.r.t massive body           [km/s]
%   params      [struct]    =   Structure containing relevant parameters
%   params.AU   [1]         =   Astronomical Unit (AU) (from DE405)                     [km]
%   params.J2   [1]         =   Gravitatonal field constant of the Earth                [-]
%   params.mu   [1]         =   Planetary consant of the planet                         [km^3/s^2]
%   params.R    [1]         =   Mean radius of the planet                               [km]
%
% OUTPUT:
%   dy      [6x1]   =   Derivative of the state vector
%   dy(1)           =   Velocity along x direction w.r.t massive body           [km/s]
%   dy(2)           =   Velocity along y direction w.r.t massive body           [km/s]
%   dy(3)           =   Velocity along z direction w.r.t massive body           [km/s]
%   dy(4)           =   Acceleration along x direction w.r.t massive body       [km/s^2]
%   dy(5)           =   Acceleration along y direction w.r.t massive body       [km/s^2]
%   dy(6)           =   Acceleration along z direction w.r.t massive body       [km/s^2]
%

% Parameters

mu = vars.mu;

omega_earth = [0;0; 7.2921e-5];
% State vector
r = y(1:3);
v = y(4:6);
v_rel = v-cross(omega_earth,r);

r_norm = norm(r);

rho = density_model(r_norm,vars);

% drag PERTURBATION EFFECT
a_drag = -1/2 * vars.A_cross *vars.Cd /vars.m *rho * norm(v_rel)^2*v_rel/norm(v_rel);

% DERIVATIVE OF THE STATE VECTOR FOR THE 2BP INCLUDING PERTURBATIONS
dy = [v; (-mu/(r_norm^3))*r + a_drag];

return