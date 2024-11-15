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


%% Mass Evaluation: 

options = optimoptions('fsolve', 'Display', 'none');
lambda0 = 0.1;
g0 = 9.80665; %m/s^2

Perc_loss =0.1; % error related to computation of parameters

% Hyp: 3 stages
% Pegasus is chosen for eps and Is of engine:

% % Specific impulse: from Pegasus
% P.Is1 = 290.2*(1-Perc_loss); % [s]
% P.Is2 = 289.4*(1-Perc_loss); % [s]
% P.Is3 = 287.4*(1-Perc_loss); % [s]

% Specific Impulse: from LauncherOne
L1.Is1 = 309*(1-Perc_loss); % [s]
L1.Is2 = 328*(1-Perc_loss); % [s]

% % Specific Impulse: from Electron
% El.Is1 = 303*(1-Perc_loss); % [s]
% El.Is2 = 333*(1-Perc_loss); % [s]

% % Structural mass index: from Pegasus
% P.eps_1 = 0.08; % [-]
% P.eps_2 = 0.09; % [-]
% P.eps_3 = 0.14; % [-]

% % Structural mass index: from LauncherOne
% L1.eps_1 = 0.06; % [-]
% L1.eps_2 = 0.06; % [-]

% Structural mass index: from Electron
El.eps_1 = 0.09; % [-]
El.eps_2 = 0.11; % [-]

IS = [L1.Is1;L1.Is2];

EPSILON = [El.eps_1;El.eps_2];

M_pay = 250;

[TR,Baseline] = Preliminary_Design(IS,EPSILON,TR.Delta_V_tot,M_pay,lambda0);

%% Mass Budget:

%TR.Mtot_Inert_Budget = TR.M.Ms1 + TR.M.Ms2 + TR.M.Ms3 + TR.M_pay; % [kg] Total inert mass budget, check that we do not go over or stay in small margin (not 1 order more)
TR.Mtot_Inert_Budget = TR.M.Ms1 + TR.M.Ms2 + TR.M_pay; % [kg] Total inert mass budget, check that we do not go over or stay in small margin (not 1 order more)

M_cables = 1.058*(TR.Length^(0.25))*sqrt(TR.M.M01); % [kg] Empirical formula slides Maggi 06, structures part 1

M_avionics = 10*(TR.M.M01)^(0.361); % [kg] Empirical formula slides Maggi 06, structures part 1

M_fairing = TR.Fair.M_fair; % [kg] calculated above 


P.Engine3.Ms3 = 126; % [kg] from Northrop Grumman Orion Series, Orion 38
P.Engine2.Ms2 = 416; % [kg] from Northrop Grumman Orion Series, Orion 50
P.Engine1.Ms1 = 1369; % [kg] from Northrop Grumman Orion Series, Orion 50S
P.Engine3.Mp3 = 770; % [kg]
P.Engine2.Mp2 = 3925; % [kg]
P.Engine3.Mp3 = 15014; % [kg]
P.Engine3.T3 = 36933.58; % [N]
P.Engine2.T2 = 131462.74; % [N]
P.Engine1.T1 = 563327.232; % [N]
P.Engine3.A_exit3 = pi*(0.52578^2)/4;
P.Engine2.A_exit2 = pi*(0.86106^2)/4;
P.Engine1.A_exit1 = pi*(1.4224^2)/4;
P.Engine1.eps_c = 10; % HYP
P.Engine3.A_t3 = (pi*(0.52578^2)/4)/10;
P.Engine2.A_t2 = (pi*(0.86106^2)/4)/10;
P.Engine1.A_t1 = (pi*(1.4224^2)/4)/10;


L1.Engine1.Ms1 = 1308; % [kg]
L1.Engine2.Ms2 = 176; % [kg]
L1.Engine1.Mp1 = 20496; % [kg]
L1.Engine2.Mp2 = 2764; % [kg]
L1.Engine1.T1 = 345162.888; % [N]
L1.Engine2.T2 = 24704.632; % [N]

TR.Engine1.T1 = L1.Engine1.T1;
TR.Engine2.T2 = L1.Engine2.T2;
%TR.Engine3.T3 = P.Engine3.T3;
TR.Engine1.A_exit1 = P.Engine1.A_exit1;
TR.Engine2.A_exit2 = P.Engine2.A_exit2;
%TR.Engine3.A_exit3 = P.Engine3.A_exit3;
TR.Engine1.A_t1 = P.Engine1.A_t1;
TR.Engine2.A_t2 = P.Engine2.A_t2;
%TR.Engine3.A_t3 = P.Engine3.A_t3;

M_struct1 = (2.55*(10^-4))* TR.Engine1.T1; % [kg] Empirical formula slides Maggi 06, structures part 1
M_struct2 = (2.55*(10^-4))* TR.Engine2.T2; % [kg] Empirical formula slides Maggi 06, structures part 1
%M_struct3 = (2.55*(10^-4))* TR.Engine3.T3; % [kg] Empirical formula slides Maggi 06, structures part 1

M_engine1 = (7.81*(10^-4))* TR.Engine1.T1 + (3.37*(10^-5))*TR.Engine1.T1*(sqrt(TR.Engine1.A_exit1/TR.Engine1.A_t1)) + 59; % [kg] Empirical formula slides Maggi 06, structures part 1
M_engine2 = (7.81*(10^-4))* TR.Engine2.T2 + (3.37*(10^-5))*TR.Engine2.T2*(sqrt(TR.Engine2.A_exit2/TR.Engine2.A_t2)) + 59; % [kg] Empirical formula slides Maggi 06, structures part 1
%M_engine3 = (7.81*(10^-4))* TR.Engine3.T3 + (3.37*(10^-5))*TR.Engine3.T3*(sqrt(TR.Engine3.A_exit3/TR.Engine3.A_t3)) + 59; % [kg] Empirical formula slides Maggi 06, structures part 1

M_parachute1 = 0.1*TR.M.Ms1; %[kg] Ele's law, anche 7%, tra 7-10%
M_parachute2 = 0.1*TR.M.Ms2; %[kg] Ele's law

%M_inert_tot_real = M_cables + M_avionics + M_fairing + TR.M_pay + M_struct1 + M_struct2 + M_struct3 +  M_engine1+  M_engine2+  M_engine3 +M_parachute1 +M_parachute2;
M_inert_tot_real = M_cables + M_avionics + M_fairing + TR.M_pay + M_struct1 + M_struct2 +  M_engine1+  M_engine2+ +M_parachute1 +M_parachute2;

Delta_inert_mass = TR.Mtot_Inert_Budget - M_inert_tot_real;

if Delta_inert_mass <=0 

    fprintf('Mass Overbudget: \n Delta_inert_mass = %d \n ',abs(Delta_inert_mass));

else

figure()
bar([TR.Mtot_Inert_Budget 0; M_inert_tot_real Delta_inert_mass],'stacked');
grid on;
xticklabels({'Team Rocket', 'Actual Masses'});
ylabel('Inert Mass');
legend('Inert Mass', 'Available Mass',Location='bestoutside');
title('Inert Mass Budget');

end

