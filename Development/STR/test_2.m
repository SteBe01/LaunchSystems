%% iterative method to verify mass estimation and to assess optimal staging

clc;
clear;
close all;

%data from other departments:
Is = [311; 311]; %[s] stages Is
dv = 8.5; %[km/s] required dv
M.pay = 250; %[kg] payload mass
OF = 2.58; %[-] Ox/Fu ratio for LOX-RP1
n = 5; %[-] load factor
diam = 1.2; %[m] external diameter
AR = 2; %aspect ratio of oblate domes [-]
loads.acc = n*9.81; %longitudinal acceleration [m/s^2]
%stage 1
M1.OF = OF;%[-] Ox/Fu ratio
M1.motor = 315; %[kg] only motor, pumps and batteries (electron - rutherford motor) %pump-fed
M1.fairing = 0; %[kg] fairing of the first stage is nonexistent
M1.rhorp1 = 807;  %[kg/m^3] density of rp1
M1.rholox = 1140; %[kg/m^3] density of lox
mat1 = 1; % 1 for Ti, 2 for Al, 3 for Steel, 4 for Carbon Fiber %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% in future versions can be optimized the material selection in function
press1 = 2; % 0 for unpressurized, 1 for pressure-fed, 2 for pump-fed, 3 for blowdown
%stage 2
M2.OF = OF;%[-] Ox/Fu ratio
M2.motor = 45; %[kg] only motor, pumps and batteries (electron - rutherford motor) %pump-fed
M2.fairing = 107; %[kg] fairing of the second stage (31.8 / 31.9)
M2.rhorp1 = 807;  %[kg/m^3] density of rp1
M2.rholox = 1140; %[kg/m^3] density of lox
mat2 = 1; % 1 for Ti, 2 for Al, 3 for Steel, 4 for Carbon Fiber %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% in future versions can be optimized the material selection in function
press2 = 2; % 0 for unpressurized, 1 for pressure-fed, 2 for pump-fed, 3 for blowdown

%first guesses:
eps0 = [0.1; 0.1]; %[0.06; 0.2]; %[-] stages structural mass indexes

%while loop parameters:
i = 2;
Nmax = 10;
err = 1;
tol = 1e-7;
eps_real = zeros(2, Nmax-1);
eps_real(:,1) = eps0;
h = 1e-5;
plotcase = 0;

while i < Nmax && err > tol
    [m_stag, m_tot, m_prop] = TANDEM(Is, eps_real(:, i-1), dv, M.pay, 0);

    %recover GLOM data:
    M1.prop = m_prop(1);
    M2.prop = m_prop(2);

    %stage 2:
    [M2, h2, th2] = inert_mass(M2, diam, AR, loads, mat2, press2);

    %stage 1:
    [M1, h1, th1] = inert_mass(M1, diam, AR, loads, mat1, press1);

    %recover eps_real:
    eps_real(:,i) = [M1.eps; M2.eps];

    %recover err:
    err = norm( eps_real(:,i) - eps_real(:,i-1) );

    %update i:
    i = i+1;
end

eps_real(:, i:end) = [];
eps_end = eps_real(:, end);

%% test

[m_stg, m_tot, m_prp] = TANDEM(Is, eps_end, dv, M.pay, 0);

%recover GLOM data:
M1.prop = m_prp(1);
M2.prop = m_prp(2);

%stage 2:
[M2test, h2test, th2test] = inert_mass(M2, diam, AR, loads, mat2, press2);

%stage 1:
[M1test, h1test, th1test] = inert_mass(M1, diam, AR, loads, mat1, press1);

%recover eps_real:
eps_test = [M1test.eps; M2test.eps];

%%

clc; clear; close all;

Is = [311; 311];
dv = 8.5;
c = Is*9.81/1000;
e = [0.05; 0.1];


lambda = 3.3:6;
z = log( c.*lambda -1 ) - log( ((e.*c).*lambda) );
size(z)
size(log( c.*lambda -1 ))
size(log( (e.*c.*lambda) ))
fun = @(lambda) dv - c'* (log( c.*lambda -1 ) - log( ((e.*c).*lambda) ) );
fun(lambda)'
%log( c*lambda -1 ./ (e.*c*lambda) )
f = -5:0.001:10;
figure(1);
plot(f, fun(f)); hold on; grid on;
plot([max(1./c), max(1./c)], [-5, 5]);

A = @(lambda) exp( dv/c(1) - c(2)*log((c(2).*lambda-1)/(e(2)*c(2)*lambda))/c(1) );
phi = @(lambda) (c(1) - A(lambda).*e(1).*c(1)).^-1;
%close all;
f3 = max(1./c)+0.1:0.001:10;
figure(2);
plot(f3, A(f3)); grid on;

phi(-5)
phi(10)

%% Functions

% function [m_stag, m_tot, m_prop] = TANDEM(Is, e, dv, m_pay, fzeroOut)
% %This function computes the optimal mass distribution and values between
% %stages for a tandem configuration.
% % INPUTS: 
% % Is : [nx1] [s] vector of the impulses of different stages
% % e  : [nx1] [1] vector of the structural mass indexes of different stages
% % dv : [1x1] [km/s] target delta_v
% % m_pay : [1x1] [kg] payload mass
% %
% % OUTPUT:
% % m_stag : [nx1] [kg] vector of the stages total masses
% % m_tot : [1x1] [kg] total initial mass
% % m_prop : [nx1] [kg] vector of the stages propellant masses
% 
% 
% g = 9.80665; %[m/s^2]
% c = Is*g/1000; %[m/s]
% 
% n = length(c);
% 
% fun = @(lambda) dv - c'*log(lambda*c-1) + log(lambda)*sum(c) + sum(c.*log(c.*e));
% 
% if fzeroOut == 0
%     options = optimset('Display','off');
% else
%     options = optimset('Display','on');
% end
% lambda0 = min( ((1-e).*c).^-1 );
% lambda0 = 0.5;
% % lambda = fsolve(fun, lambda0, options);
% lambda = fzero(fun, lambda0, options)
% 
% %condition to have m(1), m(2) > 1, therefore to have m_stag(1), m_stag(2)>0
% if lambda < lambda0
%     lambda = lambda0;
% end
% 
% m = (lambda.*c-1)./(lambda.*c.*e);
% 
% m_stag = zeros(n, 1);   %initialize
% m_stag(n) = (m(n)-1)/(1-e(n)*m(n))*m_pay;
% 
% if n == 1
%     m_stag = m_stag(1);
% else
%     i = n-1;
%     while i >= 1
%         m_stag(i) = (m(i)-1)/(1-e(i)*m(i))*(m_pay+sum(m_stag));
%         i = i-1;
%     end
% end
% 
% 
% m_prop = m_stag.*(1-e);
% m_tot = sum(m_stag) + m_pay;
% 
% end


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

fun = @(lambda) dv - c'*log( c*lambda -1 ./ ((e.*c)*lambda) );
fun = @(lambda) dv - c'* (log( c.*lambda -1 ) - log( ((e.*c).*lambda) ) );
% fun = @(lambda) dv - dot( c, log( c.*lambda -1 ./ (e.*c.*lambda) ) );
% f = 0:0.01:10;
% plot(fun(f), f)
if fsolveOut == 0
    options = optimset('Display','off');
else
    options = optimset('Display','on');
end
((1-e).*c).^-1;
lambda0 = max( ((1-e).*c).^-1 );
lambda0 = 5;
%lambda0 = max( 1./c ) + 0.05
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
% p = MEOP;


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
M.str = M.tanks + M.motor + M.fairing; %[kg] inert mass of stage (motors and tanks)
M.eps = M.str / M.tot;
end


