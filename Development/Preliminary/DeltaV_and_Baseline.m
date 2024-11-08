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

% Total: (TR struct is Team Rocket Struct)

TR.Delta_V_tot = TR.DV.Delta_V_id + TR.DV.Delta_V_g + TR.DV.Delta_V_d + TR.DV.Delta_V_s;

%Err_DV_perc = ((abs(Delta_V_base - TR.Delta_V_tot))/TR.Delta_V_tot) * 100; % [%] Error between baseline and computed DV

%% Baseline data:

%% Mass Evaluation: 

options = optimoptions('fsolve', 'Display', 'none');
lambda0 = 0.1;
g0 = 9.80665; %m/s^2

% Computation of error of our model wrt Baseline data:

P.M01 = 23130; % [kg]
[MP_model,~,~] = initialmass_opt(290.2*g0,289.4*g0,287.4*g0,0.083562229,0.095830454,0.140625,500,8165.975324,lambda0);

P.Err_M01_P = (abs(MP_model.M01 - P.M01))/P.M01; % error of our model wrt Real Pegasus data

L1.M01 = 25854; % [kg]
[ML1_model,~,~] = initialmass_opt_2stage(309*g0,328*g0,0.06,0.06,500,8462.129169,lambda0);

L1.Err_M01_L1 = (abs(ML1_model.M01c - L1.M01))/L1.M01; % error of our model wrt Real LauncherOne data

% Our Launcher:

% Hyp: 3 stages
% Pegasus is chosen for eps and Is of engine:

% Specific impulse: from Pegasus
Is1 = 290.2; % [s]
Is2 = 289.4; % [s]
Is3 = 287.4; % [s]

% % Specific Impulse: from LauncherOne
% Is1 = 309; % [s]
% Is2 = 328; % [s]

% Structural mass index: from Pegasus
eps_1 = 0.08; % [-]
eps_2 = 0.09; % [-]
eps_3 = 0.14; % [-]

% % Structural mass index: from LauncherOne
% eps_1 = 0.06; % [-]
% eps_2 = 0.06; % [-]


TR.M_pay = 400; % HELP [kg] Payload mass, do we use maximum: 400 or do we put 250??? HELP

c1 = Is1*g0; % [m/s]
c2 = Is2*g0; % [m/s]
c3 = Is3*g0; % [m/s]

% Computation of our GLOM with the assumed values:

[TR.M,TR.n,lambda_a] = initialmass_opt(c1,c2,c3,eps_1,eps_2,eps_3,TR.M_pay,TR.Delta_V_tot,lambda0);

% TR.M.M01 = TR.M.M01/(1 - P.Err_M01_P); % [kg] Corrected mass of our launcher according to erro of our model
% TR.M.M01 = TR.M.M01/(1 - L1.Err_M01_L1);  % [kg] Corrected mass of our launcher according to erro of our model

% Plot function to see min:

f = @(x) c1*log(c1*x -1) + c2*log(c2*x -1) + c3*log(c3*x -1) - log(x)*(c1+c2+c3) - c1*log(c1*eps_1) - c2*log(c2*eps_2) - c3*log(c3*eps_3) - TR.Delta_V_tot;
lambda_vec = ((10^-4):0.000001:(10^-3));
lambda_vec = lambda_vec';
y = zeros(length(lambda_vec),1);
for i = 1:length(lambda_vec)
y(i) = abs(f(lambda_vec(i)));
end

figure()
plot(lambda_vec,y);
hold on
plot(lambda_a,f(lambda_a),'o');
xlabel('\lambda [-]');
ylabel('Cost Function');
legend('Cost Function','Fsolve solution');

%% Fairing Volume Evaluation:

% Electron: payload user guide 7.0, pg 18,23; Carbon composite fairing

El.D_fair_base_int_El = 1.070; % [m]
El.D_fair_base_ext_El = 1.200; % [m]
El.H_fair_base_El = 0.566 + 0.385 + 0.0705; % [m]
El.D_fair_conetrap_int_El = 0.2783;
El.D_fair_conetrap_ext_El = 0.2783 + (0.13);
El.H_fair_conetrap_int_El = 1.3586; 
El.H_fair_nose_El = 0.119;
El.M_fair_El = 44; % [kg]
El.M_pay_max_El = 300;% [kg]

El.V_fair_El = V_fair_approx(El.D_fair_base_int_El,El.D_fair_base_ext_El,El.H_fair_base_El,El.D_fair_conetrap_int_El,El.D_fair_conetrap_ext_El,El.H_fair_conetrap_int_El,El.H_fair_nose_El);

% LauncherOne: service guide august 2020, pg 14, 15; Carbon composite fairing

L1.D_fair_base_int_L1 = 1.2624; % [m]
L1.D_fair_base_ext_L1 = 1.2624 + (2*0.11176); % [m]
L1.H_fair_base_L1 = 2.1234 + (4*0.00635); % [m]
L1.D_fair_conetrap_int_L1 = 0.44196;
L1.D_fair_conetrap_ext_L1 = 0.44196 + (2*0.11176);
L1.H_fair_conetrap_int_L1 = 3.5433-2.1234; 
L1.H_fair_nose_L1 = 3.63 - 3.5433;
%L1.M_fair_L1 = ;
L1.M_pay_max_L1 = 500; % [kg]

L1.V_fair_L1 = V_fair_approx(L1.D_fair_base_int_L1,L1.D_fair_base_ext_L1,L1.H_fair_base_L1,L1.D_fair_conetrap_int_L1,L1.D_fair_conetrap_ext_L1,L1.H_fair_conetrap_int_L1,L1.H_fair_nose_L1);

% Pegasus: Pegasus user guide 1; pg 13,15,40  
% aluminum honeycomb core with graphite/epoxy skins for the cylindrical and ogive sections of the fairing and a monocoque graphite/epoxy nose cap

P.D_fair_base_int_P = 1.153; % [m]
P.D_fair_base_ext_P = 1.27; % [m]
P.H_fair_base_P = 1.11 + 0.4; % [m]
P.D_fair_conetrap_int_P = 0.709;
P.D_fair_conetrap_ext_P = (1.27-1.153) + 0.709;
P.H_fair_conetrap_int_P = 2.656 - 1.11 + 0.4; 
P.H_fair_nose_P = 3.63 - 3.5433;
P.M_fair_P = 170; % [kg] Maggi's ex
P.M_pay_max_P = 443; % [kg]

P.V_fair_P = V_fair_approx(P.D_fair_base_int_P,P.D_fair_base_ext_P,P.H_fair_base_P,P.D_fair_conetrap_int_P,P.D_fair_conetrap_ext_P,P.H_fair_conetrap_int_P,P.H_fair_nose_P);


El.rho_pay_El = El.M_pay_max_El/El.V_fair_El.V_fair_tot;
L1.rho_pay_L1 = L1.M_pay_max_L1/L1.V_fair_L1.V_fair_tot;
P.rho_pay_P = P.M_pay_max_P/P.V_fair_P.V_fair_tot;

El.rho_pay_El_empty = El.M_pay_max_El/El.V_fair_El.V_fair_tot_int;
L1.rho_pay_L1_empty = L1.M_pay_max_L1/L1.V_fair_L1.V_fair_tot_int;
P.rho_pay_P_empty = P.M_pay_max_P/P.V_fair_P.V_fair_tot_int;

El.rho_fair_El = El.M_fair_El/(El.V_fair_El.V_fair_tot - El.V_fair_El.V_fair_tot_int);
%L1.rho_fair_L1 = L1.M_fair_L1/(L1.V_fair_L1.V_fair_tot - L1.V_fair_L1.V_fair_tot_int);
P.rho_fair_P = P.M_fair_P/(P.V_fair_P.V_fair_tot - P.V_fair_P.V_fair_tot_int);

% Volume and mass of fairing of our launcher:

TR.M_pay_max_mission = 400; % [kg]

TR.V_fair = TR.M_pay_max_mission/P.rho_pay_P; % [m^3] Whole volume
TR.V_fair_empty = TR.M_pay_max_mission/ P.rho_pay_P_empty; % [m^3] Empty volume

TR.M_fair = P.rho_fair_P * (TR.V_fair - TR.V_fair_empty); % [kg]

%% Volume of launcher: From Excel on Baseline: 

% Definition of values for baseline:

L1.fn_ratio1 = 12; % [-] fineness ratio of LauncherOne
L1.fn_ratio2 = 14; % [-] fineness ratio of LauncherOne
L1.V_empty = (0.356*(9.31+1.87)*pi*(((1.8+1.5)/2)^2)/4); % [m^3] Maggi ex
L1.V_fair_L1.Ltot = 3.63; % [m]
L1.fn_nose = L1.V_fair_L1.Ltot/L1.D_fair_base_ext_L1; % [-] Fineness ratio of nose of LauncherOne
L1.rho_launcher = 614.9798737; % [kg/m^3]
L1.V_launcher = 42.04040019; % [m^3]

P.fn_ratio = 13.31; % [-] fineness ratio of Pegasus
P.V_empty = 0.805*pi*(1.27^2)/4 + (1.905*pi*(1.27^2)/4); % [kg/m^3]
P.V_fair_P.Ltot = 2.65; % [m]
P.fn_nose = P.V_fair_P.Ltot/P.D_fair_base_ext_P; % [-] Fineness ratio of nose of Pegasus
P.rho_launcher = 1041.647642; % [kg/m^3]
P.V_launcher = 21.40839099; % [m^3]

El.V_fair_El.Ltot = 2.5; % [m]
El.fn_nose = El.V_fair_El.Ltot/El.D_fair_base_ext_El; % [-]

% Choose one baseline for the density and compute Volume of our Launcher:

TR.V_launcher = TR.M.M01/P.rho_launcher; % [m^3] Pegasus density??

TR.fn_ratio = P.fn_ratio1; % [-] L/D=f (Pegasus Baseline) imposed ??

TR.Diameter = ((4*TR.V_launcher)/(pi*TR.fn_ratio) )^(1/3);  % [m]
TR.Length = TR.fn_ratio * TR.Diameter; % [m] Whole body length

TR.fair_length = (TR.V_fair*4)/(pi*(TR.Diameter^2)); % [m] Fairing length
TR.body_length = TR.Length - TR.fair_length; % [m] Body length

% Body plot:

figure() 
plot([0 TR.fair_length TR.Length],[0 TR.Diameter TR.Diameter]/2,'k',[0 TR.fair_length TR.Length],-[0 TR.Diameter TR.Diameter]/2,'k');
axis equal
grid on
title('Geometry');
xlabel('Length [m]');
ylabel('Radius [m]');

% Plot to see where we are wrt Baseline

x_pay_eq = [El.V_fair_El.V_fair_tot;L1.V_fair_L1.V_fair_tot;P.V_fair_P.V_fair_tot];
y_pay_eq = [El.M_pay_max_El,L1.M_pay_max_L1,P.M_pay_max_P];
rho_pay_eq = polyfit(x_pay_eq,y_pay_eq,1);
y_pay_line = @(x) rho_pay_eq(1)*x + rho_pay_eq(2);
x_vec = linspace(0,10,300);

figure()
plot(x_vec,y_pay_line(x_vec));
hold on;
plot(x_pay_eq,y_pay_eq,'o');
hold on;
plot(TR.V_fair,TR.M_pay,'^');
xlabel('Volume of Fairing $[m^3]$',Interpreter='latex');
ylabel('Mass of Payload [kg]',Interpreter='latex');
title('Linear Interpolation of Payload Mass and Fairing Volume');
legend('-','Baseline','Team rocket');

x_mass_eq = [L1.V_launcher;P.V_launcher];
y_mass_eq = [L1.M01,P.M01];
V_mass_eq = polyfit(x_mass_eq,y_mass_eq,1);
y_mass_line = @(x) V_mass_eq(1)*x + V_mass_eq(2);
x_vec = linspace(0,100,300);

figure()
plot(x_vec,y_mass_line(x_vec));
hold on;
plot(x_mass_eq,y_mass_eq,'o');
hold on;
plot(TR.V_launcher,TR.M.M01,'^');
xlabel('Volume of Launcher $[m^3]$',Interpreter='latex');
ylabel('GLOM [kg]',Interpreter='latex');
title('Linear Interpolation of GLOM and Launcher Volume');
legend('-','Baseline','Team rocket');

%% Mass Budget:




%% Function for 2 stages:

function [M,n,lambda_c] = initialmass_opt_2stage(c1c,c2c,eps_1c,eps_2c,M_pay,Delta_V,lambda0)

options = optimoptions('fsolve', 'Display', 'none');


f = @(x) c1c*log(c1c*x -1) + c2c*log(c2c*x -1) - log(x)*(c1c+c2c) - c1c*log(c1c*eps_1c) - c2c*log(c2c*eps_2c) - Delta_V;

lambda_c = fsolve(f,lambda0,options);

n.n1c = (c1c*lambda_c - 1)/(c1c*eps_1c*lambda_c); % [-] 1/MR1 % must be larger than 1
n.n2c = (c2c*lambda_c - 1)/(c2c*eps_2c*lambda_c); % [-] 1/MR2

M.M2c = ((n.n2c -1)/(1 - n.n2c*eps_2c)) * (M_pay); % [kg] Mass of stage 2
M.M1c = ((n.n1c -1)/(1 - n.n1c*eps_1c)) * (M.M2c+ M_pay); % [kg] Mass of stage 1

M.Ms1c = eps_1c*M.M1c; % [kg] Mass of structure of 1
M.Ms2c = eps_2c*M.M2c; % [kg] Mass of structure of 2

M.Mp1c = M.M1c - M.Ms1c; % [kg] Mass of propellaft of 1
M.Mp2c = M.M2c - M.Ms2c; % [kg] Mass of propellant of 2

M.M02c = M_pay+ M.M2c; % [kg] Mass of stack 2
M.M01c = M.M02c + M.M1c; % [kg] Mass of stack 1


M.MR1c = 1/n.n1c;
M.MR2c = 1/n.n2c;

if n.n1c<1 | isreal(n.n1c)==0
    M.M01c = NaN;
     n.n1c = NaN;
     M.MR1c = NaN;
elseif n.n2c<1 | isreal(n.n2c)==0
 M.M01c = NaN;
  n.n2c = NaN;
  M.MR2c = NaN;

end

end
