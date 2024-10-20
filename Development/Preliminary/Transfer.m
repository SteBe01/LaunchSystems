%% Preliminary transfer

clear, clc
close all

% data

mu = astroConstants(13);

R_e = astroConstants(23);

h_airplane = 10.668e3; % [m] LauncherOne service guide august 2020

M_carrier = 0.62; % [-]  LauncherOne service guide august 2020

T_hcarr = 216.65; % [K] Wikepida temp at that h_airplane

gamma = 1.4; % [-] gamma air

R = 287; % [J/kg K] R of air

h_f = 400e3; % [m] orbit altitude by Maggi

c0 = sqrt(gamma*R*T_hcarr); % [m/s] speed of sound
v0 = M_carrier*c0; % [m/s] speed of airplane

vf = sqrt((mu*1000^3)/(h_f + (R_e* 10^3 ))); % [m/s] orbital speed

% Transfer orbit:

ra = h_f + (R_e* 10^3); % m

rp = h_airplane + (R_e* 10^3); % m

e = (ra - rp)/(ra + rp); % -
a = (ra +rp)/2; % [m]
p = a*(1-e^2); % [m]

vp = (sqrt((mu*1000^3)/p))*(1 + e); % [m/s]
va = (sqrt((mu*1000^3)/p))*(1 - e);  % [m/s]

DV1 = vp - v0;
DV2 = vf - va;

Delta_V_id = DV1 + DV2; % [m/s] Tsiolkovsky Delta V

% Losses: 

Delta_V_g = 1500; % m/s from Space mission analysis design (upper bound, medium-large missiles), also confirmed by A study of air launch methods for RLVs

Delta_V_d = 0.03 * Delta_V_id;% m/s from Space mission analysis design about 3% of total budget (upper bound)

Delta_V_s = (30.48+182.88)/2; % m/s from A study of air launch methods for RLVs btw 100 300 fps

% Total: 

Delta_V_tot = Delta_V_id + Delta_V_g + Delta_V_d + Delta_V_s;

% implement first orbit, r and v of first orbit

rr1 = [rp*10^-3 0 0]; % [km]
vv1 = [0 v0*10^-3 0]; % [km/s]
v = v0*10^-3; % [km/s]
r = rp*10^-3; % [km]

E = v^2/2 - mu/r;
a1 = -mu/(2*E);
h1 = cross(rr1,vv1);
ee = cross(vv1, h1)/mu - rr1/r;
e1 = norm(ee);

% implement final orbit, r and v of final orbit


% implement transfer orbit



% compute delta V tot


