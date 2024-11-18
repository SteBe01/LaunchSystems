
%% Tanks mass estimation
clc;
clear;
close all;

Is = [311; 311]; %[s] stages Is (electron - rutherford motor)
eps = [0.0624; 0.0304];%[0.036200000000000;0.014200000000000];%[0.058; 0.2]; %[0.035919293805387;0.014411681326543]; %[0.07; 0.12]; %[-] stages structural mass indexes
dv = 8.5; %[km/s] required dv;
M.pay = 250; %[kg] payload mass
M.fairing = 0; %[kg] fairing mass
m_motors1 = 315; %[kg] only motor, pumps and batteries (electron - rutherford motor) %pump-fed
m_motors2 = 45; %[kg] only motor, pumps and batteries (electron - rutherford motor) %pump-fed
n = 5; %[-] load factor
d = 1.2; %[m] external diameter
AR = 2; %aspect ratio of oblate domes [-]
loads.acc = n*9.81; %longitudinal acceleration [m/s^2]

%GLOM from baseline
[GLOM.m_stag, GLOM.m_tot, GLOM.m_prop] = TANDEM(Is, eps, dv, M.pay, 0);

%stage 2 analysis
M2.OF = 2.58; %[-] Ox/Fu ratio for LOX-RP1
M2.prop = GLOM.m_prop(2); %[kg] mass of propellant
% M2.rp1 = GLOM.m_prop(2) * 1 / (1+OF); %[kg] mass of rp1
% M2.lox = GLOM.m_prop(2) * OF / (1+OF);%[kg] mass of lox
% M2.prop = M2.lox + M2.rp1;%[kg] mass of propellant
M2.rhorp1 = 807;  %[kg/m^3] density of rp1
M2.rholox = 1140; %[kg/m^3] density of lox
mat2 = 1; % 1 for Ti, 2 for Al, 3 for Steel, 4 for Carbon Fiber
press2 = 2; % 0 for unpressurized, 1 for pressure-fed, 2 for pump-fed, 3 for blowdown
M2.motor = m_motors2;
M2.fairing = 0;

%find S2 info
[M2, h2, th2] = inert_mass(M2, d, AR, loads, mat2, press2);

%stage 1 analysis
M1.OF = 2.58; %[-] Ox/Fu ratio for LOX-RP1
M1.prop = GLOM.m_prop(1); %[kg] mass of propellant
% M1.rp1 = GLOM.m_prop(1) * 1 / (1+OF); %[kg] mass of rp1
% M1.lox = GLOM.m_prop(1) * OF / (1+OF);%[kg] mass of lox
% M1.prop = M1.lox + M1.rp1;%[kg] mass of propellant
M1.rhorp1 = 807;  %[kg/m^3] density of rp1
M1.rholox = 1140; %[kg/m^3] density of lox
mat1 = 1; % 1 for Ti, 2 for Al, 3 for Steel, 4 for Carbon Fiber
press1 = 2; % 0 for unpressurized, 1 for pressure-fed, 2 for pump-fed, 3 for blowdown
M1.motor = m_motors1;
M1.fairing = M.fairing;

%find S1 info
[M1, h1, th1] = inert_mass(M1, d, AR, loads, mat1, press1);

% %MA parameters
% M.M0 = M1.tot + M2.tot + m_pay + m_fairing;%[kg] initial mass (first stack mass)
% M.M0e = M.M0 - M1.prop; %[kg] mass at S1 ending
% M.MR1 = M.M0 / M.M0e; 
% M.M1 = M2.tot + m_pay + m_fairing;%[kg] second stack mass
% M.M1e = M.M1 - M2.prop; %[kg] mass at S2 ending
% M.MR2 = M.M1 / M.M1e;

eps_end = [M1.eps; M2.eps];

M.M0 = M1.tot + M2.tot;

% M.eps2 = ( M2.inert + m_pay + m_fairing ) / ( M2.tot + m_pay + m_fairing); %[-] first stage structural mass index

%% Computation of the C.O.M.:

options = optimoptions('fsolve', 'Display', 'none');
M_pay = 250;
lambda0 = 0.1;

[TR,Base] = Preliminary_Design(Is,eps,dv,M_pay,lambda0,0);

[Tank2] = tank_mass(M2, TR.Diameter, AR, loads, mat2, press2);
[Tank1] = tank_mass(M1, TR.Diameter, AR, loads, mat1, press1);

Engine2.L_nozzle=0.4; % [m]
Engine1.L_nozzle=0.4; % [m]

Engine2.m_dot_lox=5; % [kg/s]
Engine2.m_dot_fuel=5;% [kg/s]
Engine1.m_dot_lox=5;% [kg/s]
Engine1.m_dot_fuel=5;% [kg/s]

t = 500; % [s] time of flight 

[COM] = Centre_Of_Mass(TR,Tank1,Tank2,Engine1,Engine2,t);


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
lambda0 = min( ((1-e).*c).^-1 );
lambda0 = 0.5;
% lambda = fsolve(fun, lambda0, options);
lambda = fzero(fun, lambda0, options)
% lambda = real( fsolve(fun, 0.5, options) );

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

% function [M, h, t] = tank_mass(M, diam, AR, loads, mat, pressure_type)
% 
% % considers thickness equal along the whole tank.
% % safety factor to be defined.
% % the volume is the volume of propellant to be contained: you cannot use
% % this function to evaluate blowdown architectures.
% 
% %propellant masses
% mlox = M.lox; %[kg]
% mrp1 = M.rp1; %[kg]
% 
% %propellant densities
% rholox = M.rholox; %[kg/m^3]
% rhorp1 = M.rhorp1; %[kg/m^3]
% 
% %propellant volumes
% vlox = mlox / rholox; %[m^3]
% vrp1 = mrp1 / rhorp1; %[m^3]
% 
% %acceleration
% acc = loads.acc;
% 
% switch pressure_type 
%     case 0 % case for unpressurized vessel
%         MEOP = 0;
%     case 1 %case for pressure-fed
%         MEOP = 2.8 * 1e6; %[Pa] internal tank pressure
%     case 2 %case for pump-fed
%         MEOP = 700 * 1e3; %[Pa] internal tank pressure
%     case 3 % case for blowdown
%         MEOP = 50 * 1e6; %[Pa] internal tank pressure
% end
% 
% switch mat 
%     case 1 % Ti6Al4V
%         rho = 4500; %[kg/m^3]
%         t_min = 0.5 * 1e-3; %[m] minimum thickness for manufacturability
%         E = 110 * 1e9; %[Pa] young modulus
%         sy = 900 * 1e6; %[Pa] tensile yield stress
%         su = 950 * 1e6; %[Pa] tensile ultimate stress
%     case 2 % Al 2XXX
%         rho = 2700; %[kg/m^3]
%         t_min = 0.5 * 1e-3; %[m] minimum thickness for manufacturability
%         E = 70 * 1e9; %[Pa] young modulus
%         sy = 290 * 1e6; %[Pa] tensile yield stress
%         su = 390 * 1e6; %[Pa] tensile ultimate stress
%     case 3 % Steel
%         rho = 7800; %[kg/m^3]
%         t_min = 0.25 * 1e-3; %[m] minimum thickness for manufacturability
%         E = 200 * 1e9; %[Pa] young modulus
%         sy = 350 * 1e6; %[Pa] tensile yield stress
%         su = 420 * 1e6; %[Pa] tensile ultimate stress
%     case 4 % Carbon fiber
%         rho = 1800; %[kg/m^3]
%         t_min = 0.90 * 1e-3; %[m] minimum thickness for manufacturability
%         E = 250 * 1e9; %[Pa] young modulus
%         sy = 350 * 1e6; %[Pa] tensile yield stress
%         su = sy; %[Pa] tensile ultimate stress
% end
% 
% %correction factor
% Km = 1.1; 
% Kp = 1.5;
% MDP = MEOP * Km * Kp;
% jproof = 1.25;
% jburst = 1.5;
% p = MDP * jproof;
% 
% %tanks shape definition
% e = sqrt( 1 - 1 / AR^2 );         %eccentricity of oblate part [-]
% R_int = @(t) (diam - 2*t) / 2;    %internal radius of the tank [m]
% V_obl = @(t) (4/3)*pi*R_int(t)^3 / AR; %volume of the two oblate parts [m^3]
% 
% %spherical tank hypothesis
% R_sphere_lox = ( 3 * vlox / ( 4*pi ) )^(1/3); %[m] hypothetical radius of a spherical tank
% R_sphere_rp1 = ( 3 * vrp1 / ( 4*pi ) )^(1/3); %[m] hypothetical radius of a spherical tank
% 
% if R_sphere_lox < 0.98 * diam/2
%     y_lox = @(t) 2*R_sphere_lox;%fluid level inside the tank [m]
%     l_lox = @(t) y_lox(t) + 2*t;%height of the tank [m]
%     S_lox = @(t) 4 * pi * R_sphere_lox^2;%surface of tank [m^2]
% else 
%     V_cyl_lox = @(t) vlox - V_obl(t); %volume of the cylindrical part [m^3]
%     h_cyl_lox = @(t) V_cyl_lox(t) / (pi*R_int(t)^2); %height of cylindrical part [m]
%     y_lox = @(t) h_cyl_lox(t) + 2*R_int(t)/AR;%fluid level inside the tank [m]
%     l_lox = @(t) y_lox(t) + 2*t;%height of the tank [m]
%     S_cyl_lox = @(t) 2*pi*R_int(t)* h_cyl_lox(t); %surface of cylindrical part [m^2]
%     S_obl_lox = @(t) 2*pi * R_int(t)^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
%     S_lox = @(t) S_obl_lox(t) + S_cyl_lox(t); %surface of the lox tank [m^2]
% end
% 
% if R_sphere_rp1 < 0.98 * diam/2
%     y_rp1 = @(t) 2*R_sphere_rp1;%fluid level inside the tank [m]
%     l_rp1 = @(t) y_rp1(t) + 2*t;%height of the tank [m]
%     S_rp1 = @(t) 4 * pi * R_sphere_rp1^2;%surface of tank [m^2]
% else
%     V_cyl_rp1 = @(t) vrp1 - V_obl(t); %volume of the cylindrical part [m^3]
%     h_cyl_rp1 = @(t) V_cyl_rp1(t) / (pi*R_int(t)^2); %height of cylindrical part [m]
%     y_rp1 = @(t) h_cyl_rp1(t) + 2*R_int(t)/AR;%fluid level inside the tank [m]
%     l_rp1 = @(t) y_rp1(t) + 2*t;%height of the tank [m]
%     S_cyl_rp1 = @(t) 2*pi*R_int(t)* h_cyl_rp1(t); %surface of cylindrical part [m^2]
%     S_obl_rp1 = @(t) 2*pi * R_int(t)^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
%     S_rp1 = @(t) S_obl_rp1(t) + S_cyl_rp1(t); %surface of the rp1 tank [m^2]
% end
% 
% 
% %pressure at base with longitudinal acceleration (Stevin's law)
% p_lox = @(t) p + y_lox(t) * rholox * acc; %[Pa] pressure at bottom of tank during acceleration
% p_rp1 = @(t) p + y_rp1(t) * rhorp1 * acc; %[Pa] pressure at bottom of tank during acceleration
% 
% %thickness of tanks
% f_lox = @(t) t - (diam - 2*t)*p_lox(t)/( 2*sy );
% f_rp1 = @(t) t - (diam - 2*t)*p_rp1(t)/( 2*sy );
% t_lox = fzero(f_lox, 1e-3); %[m] lox tank thickness
% t_rp1 = fzero(f_rp1, 1e-3); %[m] rp1 tank thickness
% t.lox = t_lox;
% t.rp1 = t_rp1;
% 
% %check on manufacturability
% if t_lox < t_min
%     t_lox = t_min;
% elseif t_rp1 < t_min
%     t_rp1 = t_min;
% end
% 
% %height of tanks
% h.tank_lox = l_lox(t_lox); %[m] height of tank
% h.tank_rp1 = l_rp1(t_rp1); %[m] height of tank
% h.tot = h.tank_lox + h.tank_rp1; %[m] total height of tanks together
% 
% %surfaces estimation
% Slox = S_lox(t_lox);%surface of the lox tank [m^2]
% Srp1 = S_rp1(t_rp1);%surface of the rp1 tank [m^2]
% 
% %mass estimation
% M.tank_lox = rho * Slox * t_lox; %mass of the empty lox tank [kg]
% M.tank_rp1 = rho * Srp1 * t_rp1; %mass of the empty rp1 tank [kg]
% M.tot_lox = M.tank_lox + mlox; %[kg] mass of full lox tank
% M.tot_rp1 = M.tank_rp1 + mrp1; %[kg] mass of full rp1 tank
% M.tot = M.tot_lox + M.tot_rp1; %[kg] total mass of tanks and propellant
% M.tanks = M.tank_lox + M.tank_rp1; %[kg] mass of the two empty tanks
% end

function [M, h, t] = inert_mass(M, diam, AR, loads, mat, pressure_type)

% considers thickness equal along the whole tank.
% safety factor to be defined.
% the volume is the volume of propellant to be contained: you cannot use
% this function to evaluate blowdown architectures.

%propellant masses
OF = M.OF;%[-] Ox/Fu ratio
mlox = M.prop * OF / (1+OF);%[kg] mass of lox
mrp1 = M.prop * 1  / (1+OF);%[kg] mass of rp1
% mlox = M.lox; %[kg]
% mrp1 = M.rp1; %[kg]
% M.prop = mrp1 + mlox; %[kg] total propellant mass

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
        MEOP = 300 * 1e3; %[Pa] internal tank pressure
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
%p = MEOP;


%tanks shape definition
e = sqrt( 1 - 1 / AR^2 );         %eccentricity of oblate part [-]
R_int = diam / 2;    %internal radius of the tank [m]
V_obl = (4/3)*pi*R_int^3 / AR; %volume of the two oblate parts [m^3]

%spherical tank hypothesis
R_sphere_lox = ( 3 * vlox / ( 4*pi ) )^(1/3); %[m] hypothetical radius of a spherical tank
R_sphere_rp1 = ( 3 * vrp1 / ( 4*pi ) )^(1/3); %[m] hypothetical radius of a spherical tank

if R_sphere_lox < 0.98 * diam/2
    y_lox = 2*R_sphere_lox;%fluid level inside the tank [m]
    l_lox = @(t) y_lox + 2*t;%height of the tank [m]
    S_lox = 4 * pi * R_sphere_lox^2;%surface of tank [m^2]
else 
    V_cyl_lox = vlox - V_obl; %volume of the cylindrical part [m^3]
    h_cyl_lox = V_cyl_lox / (pi*R_int^2); %height of cylindrical part [m]
    y_lox = h_cyl_lox + 2*R_int/AR;%fluid level inside the tank [m]
    l_lox = @(t) y_lox + 2*t;%height of the tank [m]
    S_cyl_lox = 2*pi*R_int* h_cyl_lox; %surface of cylindrical part [m^2]
    S_obl_lox = 2*pi * R_int^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
    S_lox = S_obl_lox + S_cyl_lox; %surface of the lox tank [m^2]
end

if R_sphere_rp1 < 0.98 * diam/2
    y_rp1 = 2*R_sphere_rp1;%fluid level inside the tank [m]
    l_rp1 = @(t) y_rp1 + 2*t;%height of the tank [m]
    S_rp1 = 4 * pi * R_sphere_rp1^2;%surface of tank [m^2]
else
    V_cyl_rp1 = vrp1 - V_obl; %volume of the cylindrical part [m^3]
    h_cyl_rp1 = V_cyl_rp1 / (pi*R_int^2); %height of cylindrical part [m]
    y_rp1 = h_cyl_rp1 + 2*R_int/AR;%fluid level inside the tank [m]
    l_rp1 = @(t) y_rp1 + 2*t;%height of the tank [m]
    S_cyl_rp1 = 2*pi*R_int* h_cyl_rp1; %surface of cylindrical part [m^2]
    S_obl_rp1 = 2*pi * R_int^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
    S_rp1 = S_obl_rp1 + S_cyl_rp1; %surface of the rp1 tank [m^2]
end


%pressure at base with longitudinal acceleration (Stevin's law)
p_lox = p + y_lox * rholox * acc; %[Pa] pressure at bottom of tank during acceleration
p_rp1 = p + y_rp1 * rhorp1 * acc; %[Pa] pressure at bottom of tank during acceleration

%thickness of tanks
t_lox = diam*p_lox/( 2*sy ); %[m] lox tank thickness
t_rp1 = diam*p_rp1/( 2*sy ); %[m] rp1 tank thickness
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

%mass estimation
M.tank_lox = rho * S_lox * t_lox; %mass of the empty lox tank [kg]
M.tank_rp1 = rho * S_rp1 * t_rp1; %mass of the empty rp1 tank [kg]
M.tot_lox = M.tank_lox + mlox; %[kg] mass of full lox tank
M.tot_rp1 = M.tank_rp1 + mrp1; %[kg] mass of full rp1 tank
M.tot = M.tot_lox + M.tot_rp1 + M.motor + M.fairing; %[kg] total mass of motors, tanks and propellant
M.tanks = M.tank_lox + M.tank_rp1; %[kg] mass of the two empty tanks
M.str = M.tanks + M.motor + M.fairing;%[kg] inert mass of stage (motors and tanks)
M.eps = M.str / M.tot;
end

function [M, dv1, dv2] = ROBUST(beta, dv, dv_loss, Is, e, m_pay, h, plotcase)
%ROBUST applies the "brute force" method to compare the dv balance between
%stages to find the optimal mass for the launchers.
% First version: 09/10/2024
% Author: Pietro Bolsi
% More stages strategies can be implemented using alpha, beta vectors (dim
% = #stages -1), plot is a #stages dimensional representation.
%THEORY:
% alpha : % of ideal dv carried by the first stage
% beta : % of loss dv carried by the first stage
%
% dv_req = dv_id + dv_loss       ==dv
% dv_req = dv1 + dv2             ==dv
%
% dv1 = alpha * dv_id + beta * dv_loss
% dv2 = (1-alpha) * dv_id + (1-beta) * dv_loss
%
% Tsiolkovsky:
% dv = c * log( n ) = -c * log( MR ).  MR = Mf/M0. c = Is * g. n=1/MR.

if dv_loss > dv
    error('Losses cannot be greater than the required dv');
end

%Problem setting:
dv_id = dv - dv_loss; % [km/s]    dv == dv_req
g = 9.80665; %[m/s^2]
c = Is*g/1000; %[m/s]



if nargin < 7
    h = 0.01;
end

%Compute alpha extremes:
a_m = 1 - (1 / dv_id) * ( c(2) * log( 1/e(2) ) + ( beta-1 ) * dv_loss ) + 0.01;
a_M =     (1 / dv_id) * ( c(1) * log( 1/e(1) ) -  beta * dv_loss ) - 0.01;
if a_M < a_m
    disp('Launch not possible');
end
a_min = max(a_m, 0);
a_max = min(a_M, 1);

alpha = a_min:h:a_max ;

l = length(alpha);

%Initialize:
dv1 = zeros(l, 1); % [km/s]
dv2 = zeros(l, 1); % [km/s]
n = zeros(2, l); % [kg/kg]
M = zeros(3, l); % [kg]
M(3, :) = m_pay * ones(1, l); % [kg]
M_tot = zeros(l, 1); % [kg]
j = 0;

for i = 1:l

    a = alpha(i);

    dv1(i) = a * dv_id + beta * dv_loss;
    n(1, i) = exp( dv1(i) / c(1) );          %row containing the n of stages 1
    dv2(i) = (1-a) * dv_id + (1-beta) * dv_loss;
    n(2, i) = exp( dv2(i) / c(2) );          %row containing the n of stages 2

    M(2, i) = (n(2, i) - 1) * sum(M(:, i)) / (1 - n(2, i)*e(2)); %row containing the Mass of the Second stage wrt each alpha    
    M(1, i) = (n(1, i) - 1) * sum(M(:, i)) / (1 - n(1, i)*e(1)); %row containing the Mass of the First stage wrt each alpha
    
    M_tot(i) = sum(M(:, i));     %vector of the total masses for each alpha value

    if i > 1
        if M_tot(i) > M_tot(i-1)
            break
        end
    end
    j = j + 1;
    % M_pr2(i) = m_pay * (1-e(2)) * (1-exp(dv2(i)/c(2))) / ( exp(dv2(i)/c(2)) * e(2) - 1 );
    % M_stg2(i) = 
    % M_pr1(i) = 
end
M.tot = M_tot(j);
M_prop2 = m_pay * (1-e(2)) * (1-exp(dv2(j)/c(2))) / ( exp(dv2(j)/c(2)) * e(2) - 1 );
M_stg2 = M_prop2 / (1-e(2));
M_prop1 = M_stg2 * (1-e(1)) * (1-exp(dv1(j)/c(1))) / ( exp(dv1(j)/c(1)) * e(1) - 1 );
M_stg1 = M_prop1 / (1-e(1));
M.tot_vec = M_tot;
M.stg = [M_stg1; M_stg2];
M.prop = [M_prop1; M_prop2];


% switch plotcase
%     case 1
%     figure(1);
%     plot(alpha, M_tot);
%     xlabel('alpha');
%     ylabel('M_{tot} [kg]');
%     grid on
% 
%     figure(2);
%     plot(dv1, M_tot);
%     xlabel('Staging speed [km/s]');
%     ylabel('M_{tot} [kg]');
%     grid on
%     case 0
% end

end

