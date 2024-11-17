
%% the iterative method verifies mass estimation and assesses optimal staging

clc;
clear;
close all;

%data from other departments:
Is = [311; 311]; %[s] stages Is
dv = 8.5; %[km/s] required dv
M.pay = 250; %[kg] payload mass
OF = 2.58; %[-] Ox/Fu ratio for LOX-RP1
n = 5; %[-] load factor of longitudinal acceleration
t = 5; %[-] load factor of transversal acceleration
diam = 1.2; %[m] external diameter
AR = sqrt(3);%sqrt(3); %aspect ratio of oblate domes [-]
loads.n = n; %longitudinal acceleration [-]
loads.t = t;%transversal load factor [-]
loads.K = 2; %loads resistance safety factor [-]
Mach = 3; %[-] flight Mach number
v = 343 * Mach; %[m/s] air speed
rho_air = 0.6; %[kg/m^3] air density
maxQ = 0.5 * rho_air * v^2; %[Pa] maximum dynamic pressure
Ca = 1.3; %drag coefficient
S1 = pi * diam^2 / 4; %first stage cross section [m^2]
S2 = pi * diam^2 / 4; %second stage cross section [m^2]
loads.F_drag = maxQ * Ca * [S1-S2, S2]; %aerodynamic force acting on whole launcher [N]
h.fairing = 2; %[m] fairing height

%stage 1
M1.OF = OF;%[-] Ox/Fu ratio
M1.motor = 315; %[kg] only motor, pumps and batteries (electron - rutherford motor) %pump-fed
M1.fairing = 250; %[kg] fairing of the first stage is nonexistent
M1.rhorp1 = 807;  %[kg/m^3] density of rp1
M1.rholox = 1140; %[kg/m^3] density of lox
h1.motor = 1; %[m] height of the motor
h1.h0 = 0; %[m] starting height
mat1 = 5; % 1 for Ti, 2 for Al 2XXX, 3 for Steel, 4 for Carbon Fiber, 5 for Al 7XXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% in future versions can be optimized the material selection in function
press1 = 2; % 0 for unpressurized, 1 for pressure-fed, 2 for pump-fed, 3 for blowdown

%stage 2
M2.OF = OF;%[-] Ox/Fu ratio
M2.motor = 45; %[kg] only motor, pumps and batteries (electron - rutherford motor) %pump-fed
M2.fairing = 20; %[kg] fairing of the second stage (31.8 / 31.9)
M2.rhorp1 = 807;  %[kg/m^3] density of rp1
M2.rholox = 1140; %[kg/m^3] density of lox
h2.motor = 1; %[m] height of the motor
h2.h0 = 14; %[m] starting height
mat2 = 1; % 1 for Ti, 2 for Al 2XXX, 3 for Steel, 4 for Carbon Fiber, 5 for Al 7XXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% in future versions can be optimized the material selection in function
press2 = 2; % 0 for unpressurized, 1 for pressure-fed, 2 for pump-fed, 3 for blowdown

%first guesses:
eps0 = [0.1; 0.2]; %[0.06; 0.2]; %[-] stages structural mass indexes

%while loop parameters:
i = 2;
Nmax = 10000;
err = 1;
tol = 1e-7;
eps_real = zeros(2, Nmax-1);
eps_real(:,1) = eps0;

while i < Nmax && err > tol
    [m_stag, m_init, m_prop] = tandem_opt_staging(Is, eps_real(:, i-1), dv, M.pay, 0);

    %recover GLOM data:
    M1.prop = m_prop(1);
    M2.prop = m_prop(2);

    %stage 2:
    loads2.m = M.pay + M2.fairing;%sustained mass [kg]
    loads2.n = n;%longitudinal load factor [-]
    loads2.t = n;%transversal load factor [-]
    loads2.K = loads.K; %loads resistance safety factor [-]
    loads2.F_drag = loads.F_drag(2); %aerodynamic force [N]
    [M2, h2, th2] = inert_mass(M2, h2, diam, AR, loads2, mat2, press2);


    %stage 1:
    loads1.m = loads2.m + M2.tot;%sustained mass [kg]
    loads1.n = n;%longitudinal load factor [-]
    loads1.t = n;%transversal load factor [-]
    loads1.K = loads.K; %loads resistance safety factor [-]
    loads1.F_drag = sum( loads.F_drag ); %aerodynamic force [N]
    [M1, h1, th1] = inert_mass(M1, h1, diam, AR, loads1, mat1, press1);

    %recover eps_real:
    eps_real(:,i) = [M1.eps; M2.eps];

    %recover err:
    err = norm( eps_real(:,i) - eps_real(:,i-1) );

    %update i:
    i = i+1;
end

eps_real(:, i:end) = [];
eps_end = eps_real(:, end);

%mass related parameters
M.M0 = M.pay + M1.tot + M2.tot;%[kg] initial mass
M.M0end = M.M0 - M1.prop;
M.M1 = M.M0 - M1.tot; %[kg] mass after first stage separation
M.M1end = M.M1 - M2.prop;
M.mr1 = M.M0 / M.M0end; %[-] first stage mass ratio
M.mr2 = M.M1 / M.M1end; %[-] second stage mass ratio
M.str1 = M1.str; %[kg] first stage structural mass
M.str2 = M2.str; %[kg] second stage structural mass

%height
h.tot = h1.tot + h2.tot + h.fairing; %[m] total height
h.stg1 = h1.tot; %[m] height of first stage
h.stg2 = h2.tot; %[m] height of second stage
h.finesse_ratio = h.tot / diam; %[-] finesse ratio

%plot the two stages:
figure(1);
h2.h0 = h.stg1; %[m] updated starting height of the second stage
[~, ~, ~] = inert_mass(M2, h2, diam, AR, loads2, mat2, press2, 1);
figure(2);
[~, ~, ~] = inert_mass(M1, h1, diam, AR, loads1, mat1, press1, 1);

%% Functions:

function [m_stag, m_tot, m_prop] = tandem_opt_staging(Is, e, dv, m_pay, fzeroOut)
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

%get data
g = 9.80665; %[m/s^2]
c = Is*g/1000; %[m/s]

n = length(c);

%set root-finding problem
fun = @(k) dv - c' * log( (c-k) ./ (c.*e) );
lim = min( c-1e-2 );   %min( c.*(1-e));
if fzeroOut == 0
    options = optimset('Display','off');
else
    options = optimset('Display','on');
end
K = fzero(fun, [0, lim], options);

%get the payload ratios of each stage ( lam_i = mPL_i / mST_i )
lam = K * ( e ./ ( (1-e).*c - K ) );

%initialize
m_stag = zeros(n, 1);

if n == 1
    m_stag = m_pay / lam;
else
    i = n;
    while i >= 1
        m_stag(i) = ( m_pay + sum(m_stag) ) / lam(i) ;
        i = i-1;
    end
end
%get parameters
m_str = e .* m_stag;
m_prop = m_stag - m_str;
m_tot = sum(m_stag) + m_pay;

end

function [M, h, t] = inert_mass(M, h, diam, AR, loads, mat, pressure_type, plotcase)

% considers thickness equal along the whole tank.
% safety factor to be defined.
% the volume is the volume of propellant to be contained: you cannot use
% this function to evaluate blowdown architectures.

%constants:
g = 9.81; %[m/s^2]

%recover loads:
n = loads.n;
long_acc = n*g; %[m/s^2] longitudinal acceleration
r = loads.t;
tran_acc = r*g; %[m/s^2] transversal acceleration
m_sust = loads.m; %[kg] sustained mass

%propellant masses
OF = M.OF;%[-] Ox/Fu ratio
mlox = M.prop * OF / (1+OF);%[kg] mass of lox
mrp1 = M.prop * 1  / (1+OF);%[kg] mass of rp1

%propellant densities
rholox = M.rholox; %[kg/m^3]
rhorp1 = M.rhorp1; %[kg/m^3]

%propellant volumes
vlox = 1.10 * mlox / rholox; %[m^3] %added 10% margin
vrp1 = 1.05 * mrp1 / rhorp1; %[m^3] %added 5% margin

%recover heights:
h_motor = h.motor; %[m] height of the motor
h0 = h.h0; %[m] height of the bottom part of the stage with the launcher vertically placed on a launch pad 

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
        nu = 0.34; %[-] Poisson's ratio
    case 2 % Al 2XXX
        rho = 2700; %[kg/m^3]
        t_min = 0.5 * 1e-3; %[m] minimum thickness for manufacturability
        E = 70 * 1e9; %[Pa] young modulus
        sy = 290 * 1e6; %[Pa] tensile yield stress
        su = 390 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.33; %[-] Poisson's ratio
    case 3 % Steel
        rho = 7800; %[kg/m^3]
        t_min = 0.25 * 1e-3; %[m] minimum thickness for manufacturability
        E = 200 * 1e9; %[Pa] young modulus
        sy = 350 * 1e6; %[Pa] tensile yield stress
        su = 420 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.27; %[-] Poisson's ratio
    case 4 % Carbon fiber
        rho = 1800; %[kg/m^3]
        t_min = 0.90 * 1e-3; %[m] minimum thickness for manufacturability
        E = 250 * 1e9; %[Pa] young modulus
        sy = 350 * 1e6; %[Pa] tensile yield stress
        su = sy; %[Pa] tensile ultimate stress
        nu = 0.27; %[-] Poisson's ratio
    case 5 % Al 7XXX
        rho = 2750; %[kg/m^3]
        t_min = 1.06 * 1e-3; %[m] minimum thickness for manufacturability
        E = 70 * 1e9; %[Pa] young modulus
        sy = 500 * 1e6; %[Pa] tensile yield stress
        su = 510 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.33; %[-] Poisson's ratio
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

if R_sphere_lox < R_int
    y_lox = 2*R_sphere_lox;%fluid level inside the tank [m]
    l_lox = @(t) y_lox + 2*t;%height of the tank [m]
    S_lox = 4 * pi * R_sphere_lox^2;%surface of tank [m^2]
    h_dome_lox = R_sphere_lox; %dome height [m]
else 
    V_cyl_lox = vlox - V_obl; %volume of the cylindrical part [m^3]
    h_cyl_lox = V_cyl_lox / (pi*R_int^2); %height of cylindrical part [m]
    y_lox = h_cyl_lox + 2*R_int/AR;%fluid level inside the tank [m]
    l_lox = @(t) y_lox + 2*t;%height of the tank [m]
    S_cyl_lox = 2*pi*R_int* h_cyl_lox; %surface of cylindrical part [m^2]
    S_obl_lox = 2*pi * R_int^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
    S_lox = S_obl_lox + S_cyl_lox; %surface of the lox tank [m^2]
    h_dome_lox = R_int/AR; %dome height [m]
end

if R_sphere_rp1 < R_int
    y_rp1 = 2*R_sphere_rp1;%fluid level inside the tank [m]
    l_rp1 = @(t) y_rp1 + 2*t;%height of the tank [m]
    S_rp1 = 4 * pi * R_sphere_rp1^2;%surface of tank [m^2]
    h_dome_rp1 = R_sphere_rp1; %dome height [m]
else
    V_cyl_rp1 = vrp1 - V_obl; %volume of the cylindrical part [m^3]
    h_cyl_rp1 = V_cyl_rp1 / (pi*R_int^2); %height of cylindrical part [m]
    y_rp1 = h_cyl_rp1 + 2*R_int/AR;%fluid level inside the tank [m]
    l_rp1 = @(t) y_rp1 + 2*t;%height of the tank [m]
    S_cyl_rp1 = 2*pi*R_int* h_cyl_rp1; %surface of cylindrical part [m^2]
    S_obl_rp1 = 2*pi * R_int^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
    S_rp1 = S_obl_rp1 + S_cyl_rp1; %surface of the rp1 tank [m^2]
    h_dome_rp1 = R_int/AR; %dome height [m]
end


%pressure at base with longitudinal acceleration (Stevin's law)
p_lox = p + y_lox * rholox * long_acc; %[Pa] pressure at bottom of tank during acceleration
p_rp1 = p + y_rp1 * rhorp1 * long_acc; %[Pa] pressure at bottom of tank during acceleration

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
h.tot = h.tank_lox + h.tank_rp1 + h_motor; %[m] total height of tanks together

%plot of tanks:
if nargin > 7
    if R_sphere_lox > R_int
        bottom = @(k) h_motor + h0 + R_int/AR + t_lox -sqrt(R_int^2 - k.^2)/AR;
        top = @(k) h_motor + h0 + sqrt(R_int^2 - k.^2)/AR - R_int/AR - t_lox + h.tank_lox;
        K = linspace(-R_int, R_int, 1e4);
        plot(K, bottom(K), 'k'); grid on; axis equal; hold on;
        plot(K, top(K), 'k');
        plot([R_int, R_int], [h_motor + h0 + R_int/AR + t_lox, h_motor + h0 + R_int/AR + t_lox + h_cyl_lox], 'k');
        plot([-R_int, -R_int], [h_motor + h0 + R_int/AR + t_lox, h_motor + h0 + R_int/AR + t_lox + h_cyl_lox], 'k');
    else
        bottom = @(k) h_motor + h0 + R_sphere_lox + t_lox -sqrt(R_sphere_lox^2 - k.^2);
        top = @(k) h_motor + h0 + sqrt(R_sphere_lox^2 - k.^2) + R_sphere_lox + t_lox;
        K = linspace(-R_sphere_lox, R_sphere_lox, 1e4);
        plot(K, bottom(K), 'k'); grid on; axis equal; hold on;
        plot(K, top(K), 'k');
    end
    if R_sphere_rp1 > R_int
        bottom = @(k) h_motor + h0 + h.tank_lox + R_int/AR + t_rp1 -sqrt(R_int^2 - k.^2)/AR;
        top = @(k) h_motor + h0 + h.tank_lox + sqrt(R_int^2 - k.^2)/AR - R_int/AR - t_rp1 + h.tank_rp1;
        K = linspace(-R_int, R_int, 1e4);
        plot(K, bottom(K), 'k'); grid on; axis equal; hold on;
        plot(K, top(K), 'k');
        plot([R_int, R_int], [h_motor + h0 + R_int/AR + t_rp1 + h.tank_lox, h_motor + h0 + R_int/AR + t_rp1 + h_cyl_rp1 + h.tank_lox], 'k');
        plot([-R_int, -R_int], [h_motor + h0 + R_int/AR + t_rp1 + h.tank_lox, h_motor + h0 + R_int/AR + t_rp1 + h_cyl_rp1 + h.tank_lox], 'k');
    else
        bottom = @(k) h_motor + h0 + h.tank_lox + R_sphere_rp1 + t_rp1 -sqrt(R_sphere_rp1^2 - k.^2);
        top = @(k) h_motor + h0 + h.tank_lox + sqrt(R_sphere_rp1^2 - k.^2) + R_sphere_rp1 + t_rp1;
        K = linspace(-R_sphere_rp1, R_sphere_rp1, 1e4);
        plot(K, bottom(K), 'k'); grid on; axis equal; hold on;
        plot(K, top(K), 'k');
    end
end


%MASSES ESTIMATION

%recover material properties:
MAT.rho = rho; %[kg/m^3] material density
MAT.t_min = t_min; %[m] material manufacturability minimum thickness
MAT.E = E; %[Pa] Young modulus
MAT.sy = sy; %[Pa] yelding stress
MAT.su = su; %[Pa] ultimate stress
MAT.nu = nu; %[-] Poisson's ratio

%top connector (between top part of the stage and subsequent stage)
shape1.r = min(R_sphere_rp1, R_int);
shape1.h = h_dome_rp1 + t_rp1; 
load1.m = m_sust + M.fairing; 
load1.n = n; %longitudinal load factor
load1.K = loads.K; %safety factor
load1.F_drag = loads.F_drag; %aerodynamic drag force [N]
[Connector1.m, Connector1.th] = buckling(shape1, load1, MAT, 0);

%mass estimation of the first tank
if R_sphere_rp1 > R_int
    %validate rp1 tank size
    shape_rp1.r = R_int;
    shape_rp1.h = h_cyl_rp1;
    load_rp1.m = m_sust + M.fairing + Connector1.m + mrp1 * 0.11; %accounts for sustained masses : upper stages, first connector, fairing and rp1 tank structural mass (approximated as 11% of rp1 mass)
    load_rp1.n = n; %longitudinal load factor
    load_rp1.K = loads.K; %safety factor 
    load_rp1.F_drag = loads.F_drag; %aerodynamic drag force [N]
    [ ~, tank_rp1.th] = buckling(shape_rp1, load_rp1, MAT, MEOP);
    t_rp1 = max( t_rp1, tank_rp1.th ); %correction of rp1 tank thickness in case the previous can't sustain the load
end
M.tank_rp1 = rho * S_rp1 * t_rp1; %mass of the empty rp1 tank [kg]
M.tot_rp1 = M.tank_rp1 + mrp1; %[kg] mass of full rp1 tank

%middle connector (between two tanks)
if R_sphere_rp1 > R_int && R_sphere_lox > R_int  %that is, if tanks are cyilindrical
    shape2.r = R_int; %[m] radius of the cylindrical connector
else %that is, if at least one tank is spherical
    shape2.r = [R_sphere_rp1, R_sphere_lox, R_int];
end
shape2.h = h_dome_lox + h_dome_rp1 + t_rp1 + t_lox; %[m] lenght of the connector between the two tanks
load2.m = m_sust + M.fairing + Connector1.m + M.tot_rp1; %accounts for sustained masses : upper stages (m_sust), fairing, first connector and rp1 tank mass
load2.n = n; %longitudinal load factor
load2.K = loads.K; %safety factor
load2.F_drag = loads.F_drag; %aerodynamic drag force [N]
[Connector2.m, Connector2.th] = buckling(shape2, load2, MAT, 0);

%mass estimation of the second tank
if R_sphere_lox > R_int
    %validate lox tank size
    shape_lox.r = R_int;
    shape_lox.h = h_cyl_lox;
    load_lox.m = m_sust + M.fairing + Connector1.m + M.tot_rp1 + Connector2.m + mlox * 0.11; %accounts for sustained masses : upper stages, first and second connector, fairing, rp1 tank and lox tank structural mass (approximated as 11% of lox mass)
    load_lox.n = n; %longitudinal load factor
    load_lox.K = loads.K; %safety factor 
    load_lox.F_drag = loads.F_drag; %aerodynamic drag force [N]
    [ ~, tank_lox.th] = buckling(shape_lox, load_lox, MAT, MEOP);
    t_lox = max( t_rp1, tank_lox.th ); %correction of lox tank thickness in case the previous can't sustain the load
end
M.tank_lox = rho * S_lox * t_lox; %mass of the empty lox tank [kg]
M.tot_lox = M.tank_lox + mlox; %[kg] mass of full lox tank

%last connector (between second tank and motors)
shape3.r = min(R_sphere_lox, R_int);
shape3.h = h_dome_lox + t_lox; 
load3.m = m_sust + M.fairing + Connector1.m + M.tot_rp1 + Connector2.m + M.tot_lox; %accounts for sustained masses : upper stages, first and second connector, fairing, rp1 tank and lox tank structural mass (approximated as 11% of lox mass)
load3.n = n; %longitudinal load factor
load3.K = loads.K; %safety factor
load3.F_drag = loads.F_drag; %aerodynamic force [N]
[Connector3.m, Connector3.th] = buckling(shape3, load3, MAT, 0);

%TOTAL MASSES:
M.tot = M.tot_lox + M.tot_rp1 + M.motor + M.fairing + Connector1.m + Connector2.m + Connector3.m; %[kg] total mass of motors, tanks, propellant and connectors
M.tanks = M.tank_lox + M.tank_rp1; %[kg] mass of the two empty tanks
M.str = M.tanks + M.motor + M.fairing + Connector1.m + Connector2.m + Connector3.m; %[kg] inert mass of stage (motors, tanks, connectors and fairing)
M.eps = M.str / M.tot;
end

function [M, th] = buckling(shape, load, mat, press)

% based on NASA papers in shared folder
% computes connectors masses, heights and thicknesses to sustain
% compression loads and avoid buckling effect

%constants:
g = 9.81; %[m/s^2] gravitational acceleration

%recover loads:
m = load.m; %sustained mass [kg]
n = load.n; %longitudinal load factor [-]
K = load.K; %factor of safety [-]
F_aero = load.F_drag; %aerodynamic drag force [N]

%recover material characteristics:
E = mat.E; %[Pa] Young modulus
rho = mat.rho; %[kg/m^3]
t_min = mat.t_min; %[m] minimum thickness for manufacturability 
sy = mat.sy; %[Pa] tensile yield stress
su = mat.su; %[Pa] tensile ultimate stress
nu = mat.nu; %[-] Poisson's ratio

%recover dimensions:
r = shape.r;
h = shape.h; %distance between the base of the domes

%compute sustained load:
F_load = m * n * g + F_aero; %load [N]
if length(r) == 1 %cylindrical shell
    F_pres = pi * r^2 * press; %force from pressurant [N]
else 
    F_pres = 0;%force from pressurant [N]
end
F = F_load - F_pres; %effectively sustained load [N]

%two cases: cylindrical of trucated-cone shells:
if length(r) > 1 %trucated-cone
 
    %recover cone shape characteristics
    alpha = asin( ( r(2) - r(1) ) / h );
    L = sqrt( h^2 - (r(2) - r(1))^2 );
    l = cos(alpha) * L; %height of the shell
    r2 = cos(alpha) * r(2);
    r1 = cos(alpha) * r(1);
    
    if F < 0 %load is an axial tension load
        th = t_min;
    else %load is an axial compression load
        %compute the critical thickness for the loaded shell
        gamma = 0.33;
        th = sqrt( ( K*F*sqrt( 3*(1-nu^2) ) ) / ( 2*pi * gamma * E * cos(alpha)^2 ) ); %see page 13 of NASA paper in launch systems shared folder
    end
    %compute the mass
    S = pi * L * ( r2 + r1 ); %surface of the truncated cone
    M = S * th * rho;
else %cylindrical shell

    %recover cylinder shape characteristics
    l = h;
    
    if F < 0 %load is an axial tension load
        th = t_min;
    else %load is an axial compression load
        %compute the critical thickness for the loaded shell
        s_cr = @(t) E * ( 9 * (t/r)^1.6 + 0.16 * (t/l)^1.3 ); %[Pa] critical stress
        F_crit = @(t) 2 * pi * r * s_cr(t) * t; %[N] critical load
        f = @(t) F_crit(t) - K * F;
        th = fzero(f, [0, 1]);
    end
    %compute the mass
    M = 2*pi * r * l * th * rho; %mass of the connector
end
end

function [max_t_acc] = bending(h, loads, M, mat)

end

