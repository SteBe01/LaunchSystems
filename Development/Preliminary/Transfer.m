%% Preliminary transfer

clear, clc
close all

% data
mu = astroConstants(13);
R_e = astroConstants(23);

Delta_V_base = 8165;                            % [m/s] This one is from the Pegsus baseline to try and validate

h_airplane = 10.668e3;                          % [m] LauncherOne service guide august 2020
M_carrier = 0.62;                               % [-]  LauncherOne service guide august 2020
T_hcarr = 216.65;                               % [K] Wikepida temp at that h_airplane
gamma = 1.4;                                    % [-] gamma air
R = 287;                                        % [J/kg K] R of air
h_f = 400e3;                                    % [m] orbit altitude by Maggi

c0 = sqrt(gamma*R*T_hcarr);                     % [m/s] speed of sound
v0 = M_carrier*c0;                              % [m/s] speed of airplane

vf = sqrt((mu*1000^3)/(h_f + (R_e* 10^3 )));    % [m/s] orbital speed

% Transfer orbit:
ra = h_f + (R_e* 10^3);                         % [m]
rp = h_airplane + (R_e* 10^3);                  % [m]

e = (ra - rp)/(ra + rp);                        % [-]
a = (ra +rp)/2;                                 % [m]
p = a*(1-e^2);                                  % [m]

vp = (sqrt((mu*1000^3)/p))*(1 + e);             % [m/s]
va = (sqrt((mu*1000^3)/p))*(1 - e);             % [m/s]

DV1 = vp - v0;
DV2 = vf - va;

Delta_V_id = DV1 + DV2;                         % [m/s] Tsiolkovsky Delta V

% Losses: 
Delta_V_g = 750;                % [m/s] from Space mission analysis design (lower bound, medium-large missiles), also confirmed by A study of air launch methods for RLVs
Delta_V_d = 0.03 * Delta_V_id;  % [m/s] from Space mission analysis design about 3% of total budget (upper bound)
Delta_V_s = 30.48;              % [m/s] from A study of air launch methods for RLVs btw 100 300 fps (lower bound taken)

% Total: 
Delta_V_tot = Delta_V_id + Delta_V_g + Delta_V_d + Delta_V_s;

% Error wrt Baseline: 

Err_DV_perc = ((abs(Delta_V_base - Delta_V_tot))/Delta_V_tot) * 100; % [%] Error between baseline and computed DV

