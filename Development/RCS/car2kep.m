function [a, e, i, OM, om, th] = car2kep(r, v, mu)
% Conversion from Keplerian elements to Cartesian coordinates. Angles in radians.
%
% PROTOTYPE:
%   [a, e, i, OM, om, th] = car2kep(r, v, mu)
%
% INPUT: 
%    r       [3x1]   = position vector                [km]
%    v       [3x1]   = velocity vector                [km/s]
%    mu              = Gravitational parameter        [km^3/s^2]
% OUTPUT :
%    kepEI   [1x6]   = vector of keplerian elements
%    kepEI(1)        = a                              [km]
%    kepEI(2)        = e                              [-]
%    kepEI(3)        = i                              [rad]
%    kepEI(4)        = OM                             [rad]
%    kepEI(5)        = om                             [rad]
%    kepEI(6)        = theta                          [rad]
%
% CONTRIBUTORS:
% Camilla Airoldi, Jesus Braza Rodriguez, Elena De Marco, Riccardo Monti
%
% VERSIONS
% 2023-12-01 Last version

r_norm = norm(r); 
v_norm = norm(v);

h = cross(r, v);
h_norm = norm(h);

i = acos(h(3) / h_norm);

e_vec = 1/mu*((v_norm^2 - mu/r_norm)*r - dot(r, v)*v);
e = norm(e_vec);

eps = 1/2*v_norm^2 - mu/r_norm;
a = -mu/(2*eps);

N = cross([0 0 1], h);
N_norm = norm(N);

if(N_norm == 0)
    OM = 0;
    N = [1 0 0];
    N_norm = norm(N);
elseif(N(2) >= 0) 
    OM = acos(N(1)/N_norm);
else
    OM = 2*pi - acos(N(1)/N_norm);
end

if(e == 0)
    om = 0;
elseif(e_vec(3) >= 0) 
    om = acos(dot(N, e_vec)/(N_norm*e));
else
    om = 2*pi - acos(dot(N, e_vec)/(N_norm*e));
end

v_r = dot(r, v)/r_norm;

if(e == 0)
    th = real(acos(dot(e_vec, r)/r_norm));
elseif(v_r >= 0) 
    th = real(acos(dot(e_vec, r)/(sqrt(sum(e_vec.^2)*sum(r.^2)))));
else
    th = 2*pi - real(acos(dot(e_vec, r)/(sqrt(sum(e_vec.^2)*sum(r.^2)))));
end 

end