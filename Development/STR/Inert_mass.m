function [MASS] = Inert_mass(IS,EPSILON,Delta_V_tot,M_pay)

%% PRELIMINARY:

TR.Delta_V_tot = Delta_V_tot;
g0 = 9.80665; %m/s^2
lambda0 =0.1;

[TR,~] = Preliminary_Design(IS,EPSILON,TR.Delta_V_tot,M_pay,lambda0,0);

%% TANKS:

OF = 2.58;
MASS.OF = OF;

M1.OF = OF; %[-] Ox/Fu ratio for LOX-RP1
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
M2.OF = OF; %[-] Ox/Fu ratio for LOX-RP1
M2.prop = TR.M.Mp2; %[kg] mass of propellant
% M1.rp1 = GLOM.m_prop(1) * 1 / (1+OF); %[kg] mass of rp1
% M1.lox = GLOM.m_prop(1) * OF / (1+OF);%[kg] mass of lox
% M1.prop = M1.lox + M1.rp1;%[kg] mass of propellant
M2.rhorp1 = 807;  %[kg/m^3] density of rp1
M2.rholox = 1140; %[kg/m^3] density of lox
mat2 = 1; % 1 for Ti, 2 for Al, 3 for Steel, 4 for Carbon Fiber
press2 = 2; % 0 for unpressurized, 1 for pressure-fed, 2 for pump-fed, 3 for blowdown
M2.fairing = TR.Fair.M_fair;

AR = sqrt(3);

loads.acc = 5*g0;

[Tank2] = tank_mass(M2, TR.Diameter, AR, loads, mat2, press2);
[Tank1] = tank_mass(M1, TR.Diameter, AR, loads, mat1, press1);

%% ENGINES:

n_engines = 7;
Engine2.L_nozzle = 0.6;
Engine1.L_nozzle= 0.6;

Engine2.M = 35*7;
Engine1.M = 35*7;
Engine2.T = 27.5*(10^3)*n_engines;
Engine1.T =  27.5*(10^3)*n_engines;
n_engines = 7;

Engine2.T = 27.5*(10^3)*n_engines;
Engine1.T =  27.5*(10^3)*n_engines;

M_cables = 1.058*(TR.Length^(0.25))*sqrt(TR.M.M01); % [kg] Empirical formula slides Maggi 06, structures part 1

M_avionics = 10*(TR.M.M01)^(0.361); % [kg] Empirical formula slides Maggi 06, structures part 1

MASS.M_struct1 = (2.55*(10^-4))*Engine1.T;
MASS.M_struct2 = (2.55*(10^-4))*Engine2.T;

MASS.M_cables_distributed = M_cables/TR.Length;

Ms1_stage_ratio = TR.M.Ms1/(TR.M.Ms1+TR.M.Ms2);

MASS.M_avionics_s1 = M_avionics*Ms1_stage_ratio;

Ms2_stage_ratio = TR.M.Ms2/(TR.M.Ms1+TR.M.Ms2);

MASS.M_avionics_s2 = M_avionics*Ms2_stage_ratio;

if strcmp(Tank2.OX.FLAG, 'Cyl') && strcmp(Tank2.FU.FLAG, 'Cyl')
MASS.Diam_2 = TR.Diameter;

MASS.M_insulation_tank2_FU = 2.88*(2*pi*(Tank2.FU.H_dome*2 + Tank2.FU.H_cyl))*Tank2.FU.R_cyl_fuel;
MASS.M_insulation_tank2_OX = 1.123*(2*pi*(Tank2.OX.H_dome*2 + Tank2.OX.H_cyl))*Tank2.OX.R_cyl_lox;


elseif strcmp(Tank2.OX.FLAG, 'Sphere') && strcmp(Tank2.FU.FLAG, 'Sphere')
MASS.Diam_2 = 2*(max([Tank2.OX.R_sphere_lox;Tank2.FU.R_sphere_fuel]));

MASS.M_insulation_tank2_FU = 2.88* (4*pi*Tank2.FU.R_sphere_fuel^2);
MASS.M_insulation_tank2_OX = 1.123*(4*pi*Tank2.OX.R_sphere_lox^2);


end

if strcmp(Tank1.OX.FLAG, 'Cyl') && strcmp(Tank1.FU.FLAG, 'Cyl')
MASS.Diam_1 = TR.Diameter;

MASS.M_insulation_tank1_FU = 2.88*(2*pi*(Tank1.FU.H_dome*2 + Tank1.FU.H_cyl))*Tank1.FU.R_cyl_fuel;
MASS.M_insulation_tank1_OX = 1.123*(2*pi*(Tank1.OX.H_dome*2 + Tank1.OX.H_cyl))*Tank1.OX.R_cyl_lox;

elseif strcmp(Tank1.OX.FLAG, 'Sphere') && strcmp(Tank1.FU.FLAG, 'Sphere')
MASS.Diam_1 = 2*(max([Tank1.OX.R_sphere_lox;Tank1.FU.R_sphere_fuel]));
MASS.M_insulation_tank1_FU = 2.88* (4*pi*Tank1.FU.R_sphere_fuel^2);
MASS.M_insulation_tank1_OX = 1.123*(4*pi*Tank1.OX.R_sphere_lox^2);
end

Engine1.R_ext = MASS.Diam_1/2;
Engine2.R_ext = MASS.Diam_2/2;

MASS.Engine1 = Engine1;
MASS.Engine2 = Engine2;
MASS.Tank1 = Tank1;
MASS.Tank2 = Tank2;
MASS.TR = TR;
end % function