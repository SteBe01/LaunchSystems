function [r, v] = kep2car(kepEl, mu)
% Conversion from Keplerian elements to Cartesian coordinates. Angles in radians.
%
% PROTOTYPE:
%   [r, v] = kep2car(a, e, i, OM, om, th, mu)
%
% INPUT:
%    kepEI   [1x6]   = vector of keplerian elements
%    kepEI(1)        = a                              [km]
%    kepEI(2)        = e                              [-]
%    kepEI(3)        = i                              [rad]
%    kepEI(4)        = OM                             [rad]
%    kepEI(5)        = om                             [rad]
%    kepEI(6)        = theta                          [rad]
%   
%    mu              = Gravitational parameter        [km^3/s^2] ?????????
%
% OUTPUT:
%    r       [3x1]   = position vector                [km]
%    v       [3x1]   = velocity vector                [km/s]
%
% CONTRIBUTORS:
% Camilla Airoldi, Jesus Braza Rodriguez, Elena De Marco, Riccardo Monti
%
% VERSIONS
% 2023-12-01 Last version

[a, e, i, OM, om, th] = kepEl_demux(kepEl);

if nargin == 1
    mu = astroConstants(13);
end

p = a*(1 - e^2);
r = p/(1 + e*cos(th));

r_pf = r*[cos(th) sin(th) 0];
v_pf = sqrt(mu/p)*[-sin(th) e+cos(th) 0];

R3_OM = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
R1_i  = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)]; 
R3_om = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];

T_PF_ECI = R3_OM'*R1_i'*R3_om';

r = (T_PF_ECI*r_pf');
v = (T_PF_ECI*v_pf');

end