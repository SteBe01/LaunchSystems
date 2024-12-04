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
TR.Delta_V_tot= 10000;

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

[TR,Base] = Preliminary_Design(IS,EPSILON,TR.Delta_V_tot,M_pay,lambda0,1);

Ms2_stage_ratio = TR.M.Ms2/TR.M.Ms1;
Ms1_stage_ratio = 1- TR.M.Ms2/TR.M.Ms1;

%% Mass Budget:
%stage 1 analysis
M1.OF = 2.58; %[-] Ox/Fu ratio for LOX-RP1
M1.prop = TR.M.Mp1; %[kg] mass of propellant
% M1.rp1 = GLOM.m_prop(1) * 1 / (1+OF); %[kg] mass of rp1
% M1.lox = GLOM.m_prop(1) * OF / (1+OF);%[kg] mass of lox
% M1.prop = M1.lox + M1.rp1;%[kg] mass of propellant
M1.rhorp1 = 807;  %[kg/m^3] density of rp1
M1.rholox = 1140; %[kg/m^3] density of lox
mat1 = 1; % 1 for Ti, 2 for Al, 3 for Steel, 4 for Carbon Fiber
press1 = 2; % 0 for unpressurized, 1 for pressure-fed, 2 for pump-fed, 3 for blowdown
M1.fairing = TR.Fair.M_fair;

%stage 1 analysis
M2.OF = 2.58; %[-] Ox/Fu ratio for LOX-RP1
M2.prop = TR.M.Mp2; %[kg] mass of propellant
% M1.rp1 = GLOM.m_prop(1) * 1 / (1+OF); %[kg] mass of rp1
% M1.lox = GLOM.m_prop(1) * OF / (1+OF);%[kg] mass of lox
% M1.prop = M1.lox + M1.rp1;%[kg] mass of propellant
M2.rhorp1 = 807;  %[kg/m^3] density of rp1
M2.rholox = 1140; %[kg/m^3] density of lox
mat2 = 1; % 1 for Ti, 2 for Al, 3 for Steel, 4 for Carbon Fiber
press2 = 2; % 0 for unpressurized, 1 for pressure-fed, 2 for pump-fed, 3 for blowdown
M2.fairing = TR.Fair.M_fair;

%stage 1 analysis
M3.OF = 2.58; %[-] Ox/Fu ratio for LOX-RP1
M3.prop =TR.M.Mp1 ;%TR.M.Mp3; %[kg] mass of propellant
% M1.rp1 = GLOM.m_prop(1) * 1 / (1+OF); %[kg] mass of rp1
% M1.lox = GLOM.m_prop(1) * OF / (1+OF);%[kg] mass of lox
% M1.prop = M1.lox + M1.rp1;%[kg] mass of propellant
M3.rhorp1 = 807;  %[kg/m^3] density of rp1
M3.rholox = 1140; %[kg/m^3] density of lox
mat3 = 1; % 1 for Ti, 2 for Al, 3 for Steel, 4 for Carbon Fiber
press3 = 2; % 0 for unpressurized, 1 for pressure-fed, 2 for pump-fed, 3 for blowdown
M3.fairing = TR.Fair.M_fair;

AR = sqrt(3);
loads.acc = 5*g0;

[Tank2] = tank_mass(M2, TR.Diameter, AR, loads, mat2, press2);
[Tank1] = tank_mass(M1, TR.Diameter, AR, loads, mat1, press1);
[Tank3] = tank_mass(M3, TR.Diameter, AR, loads, mat3, press3);

t1 = 0;
t2 = 0;
t3 = 500;
tburn1 = 5000;
tburn2 = 5000;

L_nozzle2 = 0.4;
L_nozzle1 = 0.4;
L_nozzle3 = 0.4;

m_dotlox_2 = 5;
m_dotfuel_2 = 5;
m_dotlox_1 = 50;
m_dotfuel_1 = 50;
m_dotlox_3 = 2.5;
m_dotfuel_3 = 2.5;

t.t1 = t1;
t.t2 = t2;
t.t3 = t3;
t.t_burn1 = tburn1;
t.t_burn2 = tburn2;

Engine2.L_nozzle = L_nozzle2;
Engine3.L_nozzle = L_nozzle3;
Engine1.L_nozzle= L_nozzle1;

Engine2.m_dot_lox = m_dotlox_2;
Engine2.m_dot_fuel=m_dotfuel_2;
Engine1.m_dot_lox=m_dotlox_1;
Engine1.m_dot_fuel=m_dotfuel_1;
Engine3.m_dot_lox=m_dotlox_3;
Engine3.m_dot_fuel=m_dotfuel_3;

Engine2.M = 35*7;
Engine1.M = 35*7;
Engine2.t_2 = 0;
Engine1.t_1 = 0;
Engine2.T = 1.3*TR.M.M01*g0;
Engine1.T = 1.3*TR.M.M01*g0;
Engine3.T = 1.3*TR.M.M01*g0;
Engine1.R_ext = 1;
Engine2.R_ext = 1;
Engine3.R_ext = 1;

PlotFlag =0;

[COM,MoI,LOAD,TR] = Centre_Of_Mass(TR,Tank1,Tank2,Tank3,Engine1,Engine2,Engine3,t,PlotFlag);

[~,~,GEOMETRY]=COM_MoI_Stk1_dummy(mp,MASS);


%% 

loads.nx = 10*g0;
loads.nz = 3*g0;

loads.Cd = 1;
loads.Cl = 2;

v= 1.2*c0;
loads.v = v;
rho = 1.16;
loads.rho = rho;
alpha = deg2rad(20);
loads.alpha = alpha;
T = 200*10^3;
loads.T = T;

LOAD.FS=1.25;

[LOAD2_Try,TR2_Try] = Pitch_up_analysis(COM,TR,loads,LOAD,1);

[LOAD3_Try,TR2_Try] = Max_Q_analysis(COM,TR,loads,LOAD,1);

thick_2_P_up = LOAD2_Try.thick_2;
thick_2_MQ = LOAD3_Try.thick_2;
thick_2 = max([LOAD2_Try.thick_2,LOAD3_Try.thick_2]);

thick_1_P_up = LOAD2_Try.thick_4;
thick_1_MQ = LOAD3_Try.thick_4;
thick_1 = max([LOAD2_Try.thick_4,LOAD3_Try.thick_4]);

if thick_2_P_up>thick_2_MQ

Ms2_REAL = LOAD2_Try.Ms2;

elseif thick_2_P_up<thick_2_MQ

    Ms2_REAL = LOAD3_Try.Ms2;

end

if thick_1_P_up>thick_1_MQ

Ms1_REAL = LOAD2_Try.Ms1;

elseif thick_1_P_up<thick_1_MQ

    Ms1_REAL = LOAD3_Try.Ms1;

end

M01_REAL =Ms1_REAL +  Ms2_REAL + M_pay + TR.M.Mp1 +  TR.M.Mp2;

A_ref = (TR.Diameter*TR.Length)*pi;
L = (1/2)*rho*(v^2)*A_ref*loads.Cl;
D = (1/2)*rho*(v^2)*A_ref*loads.Cd;
q_D = D/TR.Length;

q_L = L/TR.Length;

b1 = LOAD.Stk1.Cone.L;

b2 = LOAD.Stk1.Stage2.L;

b3 = LOAD.Stk1.Interstage.L;

b4 = LOAD.Stk1.Stage1.L;

L_1 = q_L*b1;
L_2 = q_L*b2;
L_3 = q_L*b3;
L_4 = q_L*b4;
D_1 = q_D*b1;
D_2 = q_D*b2;
D_3 = q_D*b3;
D_4 = q_D*b4;

L_tot = L_1 + L_2 + L_3 + L_4;
D_tot = D_1 + D_2 + D_3 + D_4;

L_tot==L
D_tot==D

