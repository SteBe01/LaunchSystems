clear; clc; close all;

Delta_V_base = 8165;                            % [m/s] This one is from the Pegsus baseline to try and validate

h_carrier = 10.668e3; % m 

M_carrier = 0.62;

T_hcarr = 216.65; % K

gamma = 1.4;

R = 287; % J/kg K

h_f = 400e3; % m

c0 = sqrt(gamma*R*T_hcarr);
v0 = M_carrier*c0;

vf = sqrt((astroConstants(13)*1000^3)/(h_f + (astroConstants(23)* 10^3 ))); % m/s

% Case 1: Simple Delta V:

Delta_V_id_1 = vf - v0; % m/s

% Case 2: Hohmann transf: 

ra = h_f + (astroConstants(23)* 10^3); % m

rp = h_carrier + (astroConstants(23)* 10^3); % m

e = (ra - rp)/(ra + rp); % -
a = (ra +rp)/2; % m
p = a*(1-e^2); % m

vp = (sqrt((astroConstants(13)*1000^3)/p))*(1 + e);
va = (sqrt((astroConstants(13)*1000^3)/p))*(1 - e);

DV1 = vp - v0;
DV2 = vf - va;

Delta_V_id_2 = DV1 + DV2;

Err = Delta_V_id_1/Delta_V_id_2;

TR.DV.Delta_V_id = Delta_V_id_1; 

% Losses: 

TR.DV.Delta_V_g = 750; % m/s from Space mission analysis design, also confirmed by A study of air launch methods for RLVs

TR.DV.Delta_V_d = 0.03 * TR.DV.Delta_V_id;% m/s from Space mission analysis design about 3% of total budget

TR.DV.Delta_V_s = 30.48; % m/s from A study of air launch methods for RLVs btw 100 300 fps-> 30.48 91.44 m/s

%TR.DV.Delta_V_e = 350; % m/s speed of equator

% Total: (TR struct is Team Rocket Struct)

TR.Delta_V_tot = TR.DV.Delta_V_id + TR.DV.Delta_V_g + TR.DV.Delta_V_d + TR.DV.Delta_V_s;
Delta_V_tot= 10000;

%% Mass Evaluation: 

options = optimoptions('fsolve', 'Display', 'none');
lambda0 = 0.1;
g0 = 9.80665; %m/s^2

Perc_loss =0.1; % error related to computation of parameters

% Hyp: 3 stages
% Pegasus is chosen for eps and Is of engine:

% % Specific impulse: from Pegasus
% Base.P.Is1 = 290.2*(1-Perc_loss); % [s]
% Base.P.Is2 = 289.4*(1-Perc_loss); % [s]
% Base.P.Is3 = 287.4*(1-Perc_loss); % [s]

% Specific Impulse: from LauncherOne
Base.L1.Is1 = 309; % [s]
Base.L1.Is2 = 328; % [s]

% % Specific Impulse: from Electron
% Base.El.Is1 = 303*(1-Perc_loss); % [s]
% Base.El.Is2 = 333*(1-Perc_loss); % [s]

% % Structural mass index: from Pegasus
% Base.P.eps_1 = 0.08; % [-]
% Base.P.eps_2 = 0.09; % [-]
% Base.P.eps_3 = 0.14; % [-]

% % Structural mass index: from LauncherOne
% Base.L1.eps_1 = 0.06; % [-]
% Base.L1.eps_2 = 0.06; % [-]

% Structural mass index: from Electron
Base.El.eps_1 = 0.09; % [-]
Base.El.eps_2 = 0.11; % [-]

IS = [Base.L1.Is1;Base.L1.Is2];

EPSILON = [Base.El.eps_1;Base.El.eps_2];

M_pay = 250;

[MASS] = Inert_mass(IS,EPSILON,Delta_V_tot,M_pay);

M_prop_t1 = 1.5894e+04;
%M_prop_t1 = 0;

[X_COM1,J_y1] = COM_MoI_Stk1(M_prop_t1,MASS);

M_prop_t2 = 1.6479e+03;
%M_prop_t2 = 0;

[X_COM2,J_y2] = COM_MoI_Stk2(M_prop_t2,MASS);





