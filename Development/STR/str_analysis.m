
%% Tanks mass estimation
clc;
clear;
close all;

Is = [290; 300]; %[s] stages Is
eps = [0.07; 0.12]; %[-] stages structural mass indexes
dv = 8.5; %[m/s] required dv;
m_pay = 250; %[kg] payload mass
m_fairing = 107; %[kg] fairing mass

%GLOM from baseline
[GLOM.m_stag, GLOM.m_tot, GLOM.m_prop] = TANDEM(Is, eps, dv, m_pay, 0);

%stage 1 analysis
OF = 2.58; %[-] Ox/Fu ratio for LOX-RP1
M1.rp1 = GLOM.m_prop(1) * 1 / (1+OF); %[kg] mass of rp1
M1.lox = GLOM.m_prop(1) * OF / (1+OF);%[kg] mass of lox
M1.prop = M1.lox + M1.rp1;%[kg] mass of propellant
M1.rhorp1 = 807;  %[kg/m^3] density of rp1
M1.rholox = 1140; %[kg/m^3] density of lox

n = 10; %[-] load factor
d = 1.2; %[m] external diameter
AR = 2; %aspect ratio of oblate domes [-]
loads.acc = n*9.81; %longitudinal acceleration [m/s^2]
mat1 = 1; % 1 for Ti, 2 for Al, 3 for Steel, 4 for Carbon Fiber
press1 = 2; % 0 for unpressurized, 1 for pressure-fed, 2 for pump-fed, 3 for blowdown

%find S1 info
[M1, h1, th1] = tank_mass(M1, d, AR, loads, mat1, press1);

%stage 2 analysis
OF = 2.58; %[-] Ox/Fu ratio for LOX-RP1
M2.rp1 = GLOM.m_prop(2) * 1 / (1+OF); %[kg] mass of rp1
M2.lox = GLOM.m_prop(2) * OF / (1+OF);%[kg] mass of lox
M2.prop = M2.lox + M2.rp1;%[kg] mass of propellant
M2.rhorp1 = 807;  %[kg/m^3] density of rp1
M2.rholox = 1140; %[kg/m^3] density of lox
mat2 = 1; % 1 for Ti, 2 for Al, 3 for Steel, 4 for Carbon Fiber
press2 = 1; % 0 for unpressurized, 1 for pressure-fed, 2 for pump-fed, 3 for blowdown

%find S2 info
[M2, h2, th2] = tank_mass(M2, d, AR, loads, mat2, press2);

%MA parameters
M.M0 = M1.tot + M2.tot + m_pay + m_fairing;%[kg] initial mass (first stack mass)
M.M0e = M.M0 - M1.prop; %[kg] mass at S1 ending
M.MR1 = M.M0 / M.M0e;
M.M1 = M2.tot + m_pay + m_fairing;%[kg] second stack mass
M.M1e = M.M1 - M2.prop; %[kg] mass at S2 ending
M.MR2 = M.M1 / M.M1e;


%% Functions

function [m_stag, m_tot, m_prop] = TANDEM(Is, e, dv, m_pay, fsolveOut)
%This function computes the optimal mass distribution and values between
%stages for a tandem configuration.
% INPUTS: 
% Is : [nx1] [s] vector of the impulses of different stages
% e  : [nx1] [1] vector of the structural mass indexes of different stages
% dv : [1x1] [km/s] target delta_v
% m_pay : [1x1] [kg] payload mass
%
% OUTPUT:
% m_stag : [nx1] [kg] vector of the stages total masses
% m_tot : [1x1] [kg] total initial mass
% m_prop : [nx1] [kg] vector of the stages propellant masses


g = 9.80665; %[m/s^2]
c = Is*g/1000; %[m/s]

n = length(c);

fun = @(lambda) dv - c'*log(lambda*c-1) + log(lambda)*sum(c) + sum(c.*log(c.*e));

if fsolveOut == 0
    options = optimset('Display','off');
else
    options = optimset('Display','on');
end
lambda = fsolve(fun, 1, options);

m = (lambda.*c-1)./(lambda.*c.*e);

m_stag = zeros(n, 1);   %initialize
m_stag(n) = (m(n)-1)/(1-e(n)*m(n))*m_pay;

if n == 1

    m_stag = m_stag(1);

else
    i = n-1;
    while i >= 1

        m_stag(i) = (m(i)-1)/(1-e(i)*m(i))*(m_pay+sum(m_stag));

        i = i-1;

    end
end

m_prop = m_stag.*(1-e);
m_tot = sum(m_stag) + m_pay;

end

function [M, h, t] = tank_mass(M, diam, AR, loads, mat, pressure_type)
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
p = MDP * jproof;

%tanks shape definition
e = sqrt( 1 - 1 / AR^2 );         %eccentricity of oblate part [-]
R_int = @(t) (diam - 2*t) / 2;    %internal radius of the tank [m]
V_obl = @(t) (4/3)*pi*R_int(t)^3 / AR; %volume of the two oblate parts [m^3]

%spherical tank hypothesis
R_sphere_lox = ( 3 * vlox / ( 4*pi ) )^(1/3); %[m] hypothetical radius of a spherical tank
R_sphere_rp1 = ( 3 * vrp1 / ( 4*pi ) )^(1/3); %[m] hypothetical radius of a spherical tank

if R_sphere_lox < 0.98 * diam/2
    y_lox = @(t) 2*R_sphere_lox;%fluid level inside the tank [m]
    l_lox = @(t) y_lox(t) + 2*t;%height of the tank [m]
    S_lox = @(t) 4 * pi * R_sphere_lox^2;%surface of tank [m^2]
else 
    V_cyl_lox = @(t) vlox - V_obl(t); %volume of the cylindrical part [m^3]
    h_cyl_lox = @(t) V_cyl_lox(t) / (pi*R_int(t)^2); %height of cylindrical part [m]
    y_lox = @(t) h_cyl_lox(t) + 2*R_int(t)/AR;%fluid level inside the tank [m]
    l_lox = @(t) y_lox(t) + 2*t;%height of the tank [m]
    S_cyl_lox = @(t) 2*pi*R_int(t)* h_cyl_lox(t); %surface of cylindrical part [m^2]
    S_obl_lox = @(t) 2*pi * R_int(t)^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
    S_lox = @(t) S_obl_lox(t) + S_cyl_lox(t); %surface of the lox tank [m^2]
end

if R_sphere_rp1 < 0.98 * diam/2
    y_rp1 = @(t) 2*R_sphere_rp1;%fluid level inside the tank [m]
    l_rp1 = @(t) y_rp1(t) + 2*t;%height of the tank [m]
    S_rp1 = @(t) 4 * pi * R_sphere_rp1^2;%surface of tank [m^2]
else
    V_cyl_rp1 = @(t) vrp1 - V_obl(t); %volume of the cylindrical part [m^3]
    h_cyl_rp1 = @(t) V_cyl_rp1(t) / (pi*R_int(t)^2); %height of cylindrical part [m]
    y_rp1 = @(t) h_cyl_rp1(t) + 2*R_int(t)/AR;%fluid level inside the tank [m]
    l_rp1 = @(t) y_rp1(t) + 2*t;%height of the tank [m]
    S_cyl_rp1 = @(t) 2*pi*R_int(t)* h_cyl_rp1(t); %surface of cylindrical part [m^2]
    S_obl_rp1 = @(t) 2*pi * R_int(t)^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
    S_rp1 = @(t) S_obl_rp1(t) + S_cyl_rp1(t); %surface of the rp1 tank [m^2]
end


%pressure at base with longitudinal acceleration (Stevin's law)
p_lox = @(t) p + y_lox(t) * rholox * acc; %[Pa] pressure at bottom of tank during acceleration
p_rp1 = @(t) p + y_rp1(t) * rhorp1 * acc; %[Pa] pressure at bottom of tank during acceleration

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
h.tank_lox = l_lox(t_lox); %[m] height of tank
h.tank_rp1 = l_rp1(t_rp1); %[m] height of tank
h.tot = h.tank_lox + h.tank_rp1; %[m] total height of tanks together

%surfaces estimation
Slox = S_lox(t_lox);%surface of the lox tank [m^2]
Srp1 = S_rp1(t_rp1);%surface of the rp1 tank [m^2]

%mass estimation
M.tank_lox = rho * Slox * t_lox; %mass of the empty lox tank [kg]
M.tank_rp1 = rho * Srp1 * t_rp1; %mass of the empty rp1 tank [kg]
M.tot_lox = M.tank_lox + mlox; %[kg] mass of full lox tank
M.tot_rp1 = M.tank_rp1 + mrp1; %[kg] mass of full rp1 tank
M.tot = M.tot_lox + M.tot_rp1; %[kg] total mass of tanks and propellant
M.tanks = M.tank_lox + M.tank_rp1; %[kg] mass of the two empty tanks
end