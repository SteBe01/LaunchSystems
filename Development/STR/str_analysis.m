
%% Tanks mass estimation
clc;
clear;
close all;

M.rp1 = 2860; %[kg] mass of rp1
M.lox = 7140; %[kg] mass of lox

M.rhorp1 = 807;  %[kg/m^3] density of rp1
M.rholox = 1140; %[kg/m^3] density of lox

d = 1.2; %[m] external diameter
AR = 2; %aspect ratio of oblate domes [-]
loads.acc = 10*9.81; %longitudinal acceleration [m/s^2]
mat = 4; % 1 for Ti, 2 for Al, 3 for Steel, 4 for Carbon Fiber
press = 2; % 0 for unpressurized, 1 for pressure-fed, 2 for pump-fed, 3 for blowdown

[M, h, th] = tank_mass(M, d, AR, loads, mat, press);

%% Functions

function [M, h_tot, t] = tank_mass(M, diam, AR, loads, mat, pressure_type)
% considers thickness equal along the whole tank.
% safety factor to be defined.
% the volume is the volume of propellant to be contained: you cannot use
% this function to evaluate blowdown architectures.

%propellant masses
mlox = M.lox; %[kg]
mrp1 = M.rp1; %[kg]

%propellant densities
rholox = M.rholox; %[kg/m^3]
rhorp1 = M.rhorp1; %[kg/m^3]

%propellant volumes
vlox = mlox / rholox; %[m^3]
vrp1 = mrp1 / rhorp1; %[m^3]

%acceleration
acc = loads.acc;

switch pressure_type 
    case 0 % case for unpressurized vessel
        MEOP = 0;
    case 1 %case for pressure-fed
        MEOP = 2.8 * 1e6; %[Pa] internal tank pressure
    case 2 %case for pump-fed
        MEOP = 700 * 1e3; %[Pa] internal tank pressure
    case 3 % case for blowdown
        MEOP = 50 * 1e6; %[Pa] internal tank pressure
end

switch mat 
    case 1 % Ti6Al4V
        rho = 4500; %[kg/m^3]
        t_min = 0.5 * 1e-3; %[m] minimum thickness for manufacturability
        E = 110 * 1e9; %[Pa] young modulus
        sy = 900 * 1e6; %[Pa] tensile yield stress
        su = 950 * 1e6; %[Pa] tensile ultimate stress
    case 2 % Al 2XXX
        rho = 2700; %[kg/m^3]
        t_min = 0.5 * 1e-3; %[m] minimum thickness for manufacturability
        E = 70 * 1e9; %[Pa] young modulus
        sy = 290 * 1e6; %[Pa] tensile yield stress
        su = 390 * 1e6; %[Pa] tensile ultimate stress
    case 3 % Steel
        rho = 7800; %[kg/m^3]
        t_min = 0.25 * 1e-3; %[m] minimum thickness for manufacturability
        E = 200 * 1e9; %[Pa] young modulus
        sy = 350 * 1e6; %[Pa] tensile yield stress
        su = 420 * 1e6; %[Pa] tensile ultimate stress
    case 4 % Carbon fiber
        rho = 1800; %[kg/m^3]
        t_min = 0.90 * 1e-3; %[m] minimum thickness for manufacturability
        E = 250 * 1e9; %[Pa] young modulus
        sy = 350 * 1e6; %[Pa] tensile yield stress
        su = sy; %[Pa] tensile ultimate stress
end

%correction factor
Km = 1.1; 
Kp = 1.5;
MDP = MEOP * Km * Kp;
jproof = 1.25;
jburst = 1.5;
p = MDP * jburst;

%tanks shape definition
e = sqrt( 1 - 1 / AR^2 );         %eccentricity of oblate part [-]
R_int = @(t) (diam - 2*t) / 2;    %internal radius of the tank [m]
V_obl = @(t) (4/3)*pi*R_int(t)^3 / AR; %volume of the two oblate parts [m^3]
V_cyl_lox = @(t) vlox - V_obl(t); %volume of the cylindrical part [m^3]
V_cyl_rp1 = @(t) vrp1 - V_obl(t); %volume of the cylindrical part [m^3]
h_cyl_lox = @(t) V_cyl_lox(t) / (pi*R_int(t)^2); %height of cylindrical part [m]
h_cyl_rp1 = @(t) V_cyl_rp1(t) / (pi*R_int(t)^2); %height of cylindrical part [m]
l_lox = @(t) h_cyl_lox(t) + 2*R_int(t)/AR + 2*t; %height of the tank [m]
l_rp1 = @(t) h_cyl_rp1(t) + 2*R_int(t)/AR + 2*t; %height of the tank [m]

%pressure at base with longitudinal acceleration (Stevin's law)
p_lox = @(t) p + l_lox(t) * rholox * acc; %[Pa] pressure at bottom of tank during acceleration
p_rp1 = @(t) p + l_rp1(t) * rhorp1 * acc; %[Pa] pressure at bottom of tank during acceleration

%thickness of tanks
f_lox = @(t) t - (diam - 2*t)*p_lox(t)/( 2*sy );
f_rp1 = @(t) t - (diam - 2*t)*p_rp1(t)/( 2*sy );
t_lox = fzero(f_lox, 1e-3); %[m] lox tank thickness
t_rp1 = fzero(f_rp1, 1e-3); %[m] rp1 tank thickness
t.lox = t_lox;
t.rp1 = t_rp1;

%check on manufacturability
if t_lox < t_min
    t_lox = t_min;
elseif t_rp1 < t_min
    t_rp1 = t_min;
end

%height of tanks
h_lox = l_lox(t_lox); %[m] height of tank
h_rp1 = l_rp1(t_rp1); %[m] height of tank
h_tot = h_lox + h_rp1; %[m] total height of tanks together

%surfaces estimation
S_cyl_lox = 2*pi*R_int(t_lox)* h_cyl_lox(t_lox); %surface of cylindrical part [m^2]
S_cyl_rp1 = 2*pi*R_int(t_rp1)* h_cyl_rp1(t_rp1); %surface of cylindrical part [m^2]
S_obl_lox = 2*pi * R_int(t_lox)^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
S_obl_rp1 = 2*pi * R_int(t_rp1)^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
S_lox = S_obl_lox + S_cyl_lox; %surface of the lox tank [m^2]
S_rp1 = S_obl_rp1 + S_cyl_rp1; %surface of the rp1 tank [m^2]

%mass estimation
M.tank_lox = rho * S_lox * t_lox; %mass of the empty lox tank [kg]
M.tank_rp1 = rho * S_rp1 * t_rp1; %mass of the empty rp1 tank [kg]
M.tot_lox = M.tank_lox + mlox; %[kg] mass of full lox tank
M.tot_rp1 = M.tank_rp1 + mrp1; %[kg] mass of full rp1 tank
M.tot = M.tot_lox + M.tot_rp1; %[kg] total mass of tanks and propellant
end