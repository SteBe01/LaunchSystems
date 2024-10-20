%% Preliminary transfer

clear, clc
close all

% data

mu = 398600;

h_airplane = 10;            % [km]
v_airplane = 800/3600;      % [km/s]
h_final_orbit = 250;        % [km]
Earth_radius = 6371;        % [km]

r_airplane = Earth_radius + h_airplane;
r_final_orbit = Earth_radius + h_final_orbit;

% implement first orbit, r and v of first orbit

rr = [r_airplane 0 0];
vv = [0 v_airplane 0];
v = v_airplane;
r = r_airplane;

E = v^2/2 - mu/r;
a1 = -mu/(2*E);
h1 = cross(rr,vv);
ee = cross(vv, h1)/mu - rr/r;
e1 = norm(ee);

% implement final orbit, r and v of final orbit



% implement transfer orbit



% compute delta V tot


