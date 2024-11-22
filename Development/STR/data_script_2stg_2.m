
%% the iterative method verifies mass estimation and assesses optimal staging

clc;
clear;
close all;

%data from other departments:
Is = [328; 343]; %[s] stages Is
dv = 8.5; %[km/s] required dv
M.pay = 250; %[kg] nominal payload mass
M.pay_max = 400; %[kg] maximum payload mass
OF = 2.58; %[-] Ox/Fu ratio for LOX-RP1
nx = 5; %[-] load factor of longitudinal acceleration
ny = 5; %[-] load factor of transversal acceleration
diam = 1.1; %[m] external diameter
AR = sqrt(3);%sqrt(3); %aspect ratio of oblate domes [-]
loads.n = nx; %longitudinal acceleration [-]
loads.t = ny;%transversal load factor [-]
loads.K = 3; %loads resistance safety factor [-]
maxQ = 45000; %0.5 * rho_air * v^2; %[Pa] maximum dynamic pressure
Ca1 = 1.3; %drag coefficient of first stage
Ca2 = 1.3; %drag coefficient of second stage
Caf = 1.3; %drag coefficient of the fairing
S1 = pi * diam^2 / 4; %first stage cross section [m^2]
S2 = pi * diam^2 / 4; %second stage cross section [m^2]
loads.F_drag = zeros(3, 1); %aerodynamic forces acting on the three part of the launcher (fairing, stg2, stg1) [N]
TW1 = 1.3; %[-] T/W ratio of first stage
TW2 = 1.1; %[-] T/W ratio of second stage

%stage 1
M1.OF = OF;%[-] Ox/Fu ratio
M1.motor = 315; %[kg] only motor, pumps and batteries (electron - rutherford motor) %pump-fed
M1.rhorp1 = 807;  %[kg/m^3] density of rp1
M1.rholox = 1140; %[kg/m^3] density of lox
h1.motor = 0.89; %[m] height of the motor
h1.h0 = 0; %[m] starting height
mat1 = 2; % 1 for Ti, 2 for Al 2XXX, 3 for Steel, 4 for Carbon Fiber, 5 for Al 7XXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% in future versions can be optimized the material selection in function
press1 = 2; % 0 for unpressurized, 1 for pressure-fed, 2 for pump-fed, 3 for blowdown

%stage 2
M2.OF = OF;%[-] Ox/Fu ratio
M2.motor = 45; %[kg] only motor, pumps and batteries (electron - rutherford motor) %pump-fed
M2.rhorp1 = 807;  %[kg/m^3] density of rp1
M2.rholox = 1140; %[kg/m^3] density of lox
h2.motor = 0.89; %[m] height of the motor
h2.h0 = 14; %[m] starting height
mat2 = 2; % 1 for Ti, 2 for Al 2XXX, 3 for Steel, 4 for Carbon Fiber, 5 for Al 7XXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% in future versions can be optimized the material selection in function
press2 = 2; % 0 for unpressurized, 1 for pressure-fed, 2 for pump-fed, 3 for blowdown

%first guesses:
eps0 = [0.1; 0.2]; %[0.06; 0.2]; %[-] stages structural mass indexes
fairing.base_diam = diam; %[m] first guess for fairing base diameter
adapter.base_diam = fairing.base_diam; %[m] first guess for fairing and adapter base diameter
M.M0 = 10000; %[kg]
M.M1 = 1300; %[kg]
h.tot = 16; %[m]


%while loop parameters:
i = 2;
Nmax = 10000;
err = 1;
tol = 1e-7;
eps_real = zeros(2, Nmax-1);
eps_real(:,1) = eps0;
diam_fairing = zeros(1, Nmax-1);
diam_fairing(1) = fairing.base_diam;

while i < Nmax && err > tol
    %recover MER data:
    M.avionics = 10 * M.M0^0.361; %[kg]
    M.wiring = 1.058 * h.tot^0.25 * sqrt(M.M0); %[kg]

    %divide MER data between stages:
    M1.avionics = 0.5 * M.avionics;%[kg]
    M2.avionics = 0.5 * M.avionics;%[kg]
    M1.wiring = M.wiring * m_stag(1) / sum(m_stag);%[kg]
    M2.wiring = M.wiring * m_stag(2) / sum(m_stag);%[kg]
    M1.T_struct = 2.55 * M.M0 * 9.81 * TW1 * 1e-4; %[kg]
    M2.T_struct = 2.55 * M.M1 * 9.81 * TW2 * 1e-4; %[kg] 

    %optimal staging:
    [m_stag, m_init, m_prop] = tandem_opt_staging(Is, eps_real(:, i-1), dv, M.pay, 0);

    %recover GLOM data:
    M1.prop = m_prop(1);
    M2.prop = m_prop(2);

    %compute bending loads:
    loads.M_exp = loads.ny * M.M0 * 9.81 * h.tot / 4;

    %adapter 
    loads_a = loads; %recover loads
    loads_a.m = M.pay;%sustained mass [kg]
    loads_a.F_drag = 0; %aerodynamic load is null for the payload adapter [N]
    [adapter] = adapter_fun(adapter, loads_a);

    %fairing:
    Sf = pi * diam_fairing(i-1)^2 / 4; %fairing cross section [m^2]
    loads_f = loads; %recover loads
    loads.F_drag(3) = maxQ * Caf * Sf; %compose aerodynamic forces vector [N]
    loads_f.F_drag = loads.F_drag(3); %aerodynamic force acting on fairing [N]
    [fairing] = fairing_fun(M.pay_max, M.pay, fairing, loads_f);

    %stage 2:
    M2.R_next = diam_fairing(i-1) / 2; %[m]
    M2.fairing = adapter.m + fairing.m + M2.avionics + M2.wiring; %[kg] WE ASSUME THE ADAPTER DETATCH FROM THE LAUNCHER TOGETHER WITH THE SECOND STAGE
    loads2 = loads; %recover loads
    loads2.m = M.pay + M2.fairing;%sustained mass [kg]
    if (S2-Sf) > 0
        loads.F_drag(2) = maxQ * Ca2 * (S2-Sf); %compose aerodynamic forces vector [N]
    else
        loads.F_drag(2) = 0; %compose aerodynamic forces vector [N]
    end
    loads2.F_drag = sum( loads.F_drag ); %aerodynamic force [N]
    [M2, h2, th2] = inert_mass_common_dome(M2, h2, diam, AR, loads2, mat2, press2);

    %stage 1:
    M1.R_next = M2.R_end; %[m]
    M1.fairing = 0; %[kg] WE ASSUME THE FAIRING DETATCH FROM THE LAUNCHER TOGETHER WITH THE FIRST STAGE
    loads1 = loads; %recover loads
    loads1.m = M1.avionics + M1.wiring + M2.tot + M.pay;%sustained mass [kg] M2.tot comprises the adapter
    if (S1-S2-Sf) > 0
        loads.F_drag(1) = maxQ * Ca1 * (S1-S2-Sf); %compose aerodynamic forces vector [N]
    else
        loads.F_drag(1) = 0; %compose aerodynamic forces vector [N]
    end
    loads1.F_drag = sum( loads.F_drag ); %aerodynamic force [N]
    [M1, h1, th1] = inert_mass_common_dome(M1, h1, diam, AR, loads1, mat1, press1);

    %recover eps_real:
    eps_real(:,i) = ( eps_real(:,i-1) + [M1.eps; M2.eps] ) / 2;

    %recover real dimensions of M1, M2, of the fairing and of the adapter:
    diam_fairing(i) = M2.diam1; %[m]
    fairing.base_diam = diam_fairing(i); %[m]
    adapter.base_diam = diam_fairing(i); %[m]
    h.tot = h1.tot + h1.C1 + h1.motor + h2.tot + fairing.L; %[m] total height
    h.CG = ( h1.CG.tot * (M1.tot - M1.fairing) + h2.CG.tot * (M2.tot - M2.fairing) +...
            + ( fairing.h0 + (1/3) * fairing.L ) * (fairing.m + adapter.m + M.pay) ) / M.M0;
    M.M1 = M.M0 - M1.tot; %[kg] mass after first stage separation
    
    %recover err:
    err = norm( [eps_real(:,i); diam_fairing(i)] - [eps_real(:,i-1); diam_fairing(i-1)] );% + 1e-7 * norm( diam_fairing(i) - diam_fairing(i-1) );

    %update i:
    i = i+1;
end

eps_real(:, i:end) = [];
eps_end = eps_real(:, end);
diam_fairing(:, i:end) = [];

%recover loads
loads.F_drag_tot = sum( loads.F_drag ); %total aerodynamic drag [N]

%mass related parameters
M.M0 = M.pay + M1.tot + M2.tot;%[kg] initial mass
M.M0end = M.M0 - M1.prop;
M.M1 = M.M0 - M1.tot; %[kg] mass after first stage separation
M.M1end = M.M1 - M2.prop;
M.mr1 = M.M0 / M.M0end; %[-] first stage mass ratio
M.mr2 = M.M1 / M.M1end; %[-] second stage mass ratio
M.str1 = M1.str; %[kg] first stage structural mass
M.str2 = M2.str; %[kg] second stage structural mass
M.fairing = fairing.m; %[kg] fairing mass
M.adapter = adapter.m; %[kg] payload adapter mass

%height
h.tot = h1.tot + h1.C1 + h1.motor + h2.tot + fairing.L; %[m] total height
h.stg1 = h1.tot + h1.C1 + h1.motor; %[m] height of first stage
h.stg2 = h2.tot + h2.C1 + h2.motor; %[m] height of second stage
h.fairing = fairing.L; %[m] height of the fairing
h.finesse_ratio = h.tot / diam; %[-] finesse ratio

%plot the two stages:
figure(1);
fairing.h0 = h.tot - h.fairing;%[m] updated starting height of the fairing
[~] = fairing_fun(M.pay_max, M.pay_max, fairing, loads_f, 1);
figure(2);
h2.h0 = h.stg1 - h2.C3 - h2.motor; %[m] updated starting height of the second stage
[~, ~, ~] = inert_mass(M2, h2, diam, AR, loads2, mat2, press2, 1);
figure(3);
[~, ~, ~] = inert_mass(M1, h1, diam, AR, loads1, mat1, press1, 1);

%% Functions

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

function [M, h, th] = inert_mass_common_dome(M, h, diam, AR, loads, mat_id, pressure_type, plotcase)

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

%recover dimensions:
h_motor = h.motor; %[m] height of the motor
h0 = h.h0; %[m] height of the bottom part of the stage with the launcher vertically placed on a launch pad 
R_next = M.R_next; %[m] radius of the subsequent stage

%get material properties
mat = mat_switch(mat_id);
rho = mat.rho;%[kg/m^3] material density
t_min = mat.t_min; %[m] material manufacturability minimum thickness
E = mat.E; %[Pa] Young modulus
sy = mat.sy; %[Pa] yelding stress
su = mat.su; %[Pa] ultimate stress
nu = mat.nu; %[-] Poisson's ratio

%get MEOP (maximum expected operating pressure)
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

%correction factor
Km = 1.1; 
Kp = 2;
MDP = MEOP * Km * Kp;
jproof = 1.25;
jburst = 1.5;
p = MDP * jburst;

%tanks shape definition
e = sqrt( 1 - 1 / AR^2 );         %eccentricity of oblate part [-]
R_int = diam / 2;    %internal radius of the tank [m]
V_obl = (4/3)*pi*R_int^3 / AR; %volume of the two oblate parts [m^3]

%check on tank shape
v_min = min(vlox, vrp1);

if v_min < (4/3) * pi * R_int^3 / AR  %adopt spherical tanks
    %lox
    R_lox = ( 3 * vlox / ( 4*pi ) )^(1/3); %[m] radius of lox spherical tank
    h_dome_lox = R_lox; %dome height [m]
    h_cyl_lox = 0;%height of cylindrical part [m]
    S_lox = 4 * pi * R_lox^2;%surface of tank [m^2]
    AR_lox = 1;
    h0_lox = (2/3)*R_lox + h0 + h_motor;%[m]

    %rp1
    R_rp1 = ( 3 * vrp1 / ( 4*pi ) )^(1/3); %[m] radius of rp1 spherical tank
    h_dome_rp1 = R_rp1; %dome height [m]
    h_cyl_rp1 = 0;%height of cylindrical part [m]
    S_rp1 = 4 * pi * R_rp1^2;%surface of tank [m^2]
    y_lox = 2*R_lox;%fluid level inside the tank [m]
    AR_rp1 = 1;
    h0_rp1 = h0_lox + 2*h_dome_lox + (1/2)*R_rp1;%[m]
    h0_C1 = h0_rp1 + R_rp1;%[m]
    h0_C2 = h0_lox + R_lox;%[m]
    h0_C3 = h0 + h_motor;%[m]
else  %adopt cyl tanks with shared dome
    %lox
    R_lox = R_int; %[m] radius of the lox tank
    h_dome_lox = R_lox/AR; %dome height [m]
    V_cyl_lox = vlox - V_obl; %volume of the cylindrical part [m^3]
    h_cyl_lox = V_cyl_lox / (pi*R_lox^2); %height of cylindrical part [m]
    y_lox = h_cyl_lox + 2*R_lox/AR;%fluid level inside the tank [m]
    S_cyl_lox = 2*pi*R_lox* h_cyl_lox; %surface of cylindrical part [m^2]
    S_obl_lox = 2*pi * R_lox^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
    S_lox = S_obl_lox + S_cyl_lox; %surface of the lox tank [m^2]
    AR_lox = AR;
    h0_lox = (2/3)*R_lox + h0 + h_motor;%[m] 

    %rp1
    R_rp1 = R_int; %[m] radius of the rp1 tank
    h_dome_rp1 = R_rp1/AR; %dome height [m]
    V_cyl_rp1 = vrp1; %volume of the cylindrical part [m^3]
    h_cyl_rp1 = V_cyl_rp1 / (pi*R_rp1^2); %height of cylindrical part [m]
    y_rp1 = h_cyl_rp1 + 2*R_rp1/AR;%fluid level inside the tank [m]
    S_cyl_rp1 = 2*pi*R_rp1* h_cyl_rp1; %surface of cylindrical part [m^2]
    S_obl_rp1 = 2*pi * R_rp1^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
    S_rp1 = S_obl_rp1 + S_cyl_rp1; %surface of the rp1 tank [m^2]
    AR_rp1 = AR;
    h0_rp1 = h0_lox + h_dome_lox + h_cyl_lox + 0.16; %[m]. 0.16m added based on Edberg-Costa to account for insulation
    h0_C1 = h0_rp1 + h_cyl_rp1;%[m]
    h0_C3 = h0 + h_motor;%[m]
end

%recover dimensions
l_lox = @(t) y_lox + 2*t;%height of the tank [m]
l_rp1 = @(t) y_rp1 + 2*t;%height of the tank [m]
h.CG.lox = h_cyl_lox / 2 + h_dome_lox + h0_lox; %[m] COG of lox tank
h.CG.rp1 = h_cyl_rp1 / 2 + h_dome_rp1 + h0_rp1; %[m] COG of rp1 tank
M.R_end = R_lox / 2; %[m] ending radius of the stage

%pressure at base with longitudinal acceleration (Stevin's law)
p_lox = p + y_lox * rholox * long_acc; %[Pa] pressure at bottom of tank during acceleration
p_rp1 = p + y_rp1 * rhorp1 * long_acc; %[Pa] pressure at bottom of tank during acceleration

%thickness of tanks for pressure
t_p_lox = diam*p_lox/( 2*sy ); %[m] lox tank thickness
t_p_rp1 = diam*p_rp1/( 2*sy ); %[m] rp1 tank thickness


%MASSES ESTIMATION

%top interstage (first connector) (between top part of the stage and subsequent stage)
if R_next > R_rp1  %that is, if  next stage radius is NOT less than this one's
    shape1.r = R_rp1; %[m] radius of the cylindrical connector
    shape1.h = shape1.r * 5 / 2;
else %that is, if  next stage radius is less than this one's
    shape1.r = [R_next, R_rp1]; %for simplicity we take the same dimensions of the "both-spherical" case
    shape1.h = shape1.r(2) * 5 / 2;
end 
load1 = loads;
load1.p = 0; %internal pressure [Pa]
if nargin < 8
    [C1.m, C1.th] = buckling_bending(shape1, load1, mat);
else
    shape1.h0 = h0_C1;
    [C1.m, C1.th, C1.XY] = buckling(shape1, load1, mat, 1);
    plot(C1.XY(1,:), C1.XY(2,:), '--k', DisplayName='true');grid on; axis equal; hold on;
end
h.CG.C1 = h0_C1 + shape1.h / 2;%[m] exact solution for cyl, approx for cones


%mass estimation of the first tank
if h_cyl_rp1 > 0    %cyl tank
    %validate rp1 tank size
    shape_rp1.r = R_rp1;
    shape_rp1.h = h_cyl_rp1;
    load_rp1 = loads;
    load_rp1.m = m_sust + C1.m + mrp1 * 0.11; %accounts for sustained masses : upper stages, first connector, fairing and rp1 tank structural mass (approximated as 11% of rp1 mass)
    load_rp1.F_drag = loads.F_drag; %aerodynamic drag force [N]
    load_rp1.p = MEOP;%[Pa]
    [~, t_bb_rp1] = buckling_bending(shape_rp1, load_rp1, mat);
    t_rp1 = max( t_p_rp1, t_bb_rp1 ); %correction of rp1 tank thickness in case the previous can't sustain the load
else %spherical tank
    t_rp1 = max( t_p_rp1, t_min );
end
th.rp1 = t_rp1;%[m]
M.tank_rp1 = rho * S_rp1 * t_rp1; %mass of the empty rp1 tank [kg]
M.tot_rp1 = M.tank_rp1 + mrp1; %[kg] mass of full rp1 tank

%middle connector (intertank) (between two tanks)
if R_rp1 == R_lox %that is, if tanks have the same diameter
    C2.m = 0;%[kg]
    shape2.h = 0; %[m]
    h.CG.C2 = 0; %[m]
else %that is, if at least one tank is spherical
    shape2.r = [R_rp1, R_lox]; %for simplicity we take the same dimensions of the "both-spherical" case
    shape2.h = 0.5 * R_rp1 + h_dome_lox + h_dome_rp1;%[m] lenght of the connector between the two tanks
    load2 = loads;
    load2.m = m_sust + C1.m + M.tot_rp1; %accounts for sustained masses : upper stages (m_sust), fairing, first connector and rp1 tank mass
    load2.F_drag = loads.F_drag; %aerodynamic drag force [N]
    if nargin < 8
        [C2.m, C2.th] = buckling_bending(shape2, load2, mat);
    else
        shape2.h0 = h0_C2;
        [C2.m, C2.th, C2.XY] = buckling_bending(shape2, load2, mat, 1);
        plot(C2.XY(1,:), C2.XY(2,:), '--k', DisplayName='true');
    end
    h.CG.C2 = h0_C2 + shape2.h / 2;%[m] exact solution for cyl, approx for cones
end

%mass estimation of the second tank
if h_cyl_lox > 0 %cyl tank
    %validate lox tank size
    shape_lox.r = R_int;
    shape_lox.h = h_cyl_lox;
    load_lox = loads;
    load_lox.m = m_sust + C1.m + M.tot_rp1 + C2.m + mlox * 0.11; %accounts for sustained masses : upper stages, first and second connector, fairing, rp1 tank and lox tank structural mass (approximated as 11% of lox mass)
    load_lox.F_drag = loads.F_drag; %aerodynamic drag force [N]
    [t_bb_lox] = buckling_bending(shape_lox, load_lox, mat, MEOP);
    t_lox = max( t_p_lox, t_bb_lox ); %correction of lox tank thickness in case the previous can't sustain the load
else %spherical tank
    t_lox = max( t_p_lox, t_min);
end
th.lox = t_lox;%[m]
M.lox_insulation = 1.12 * S_lox; %[kg] mass of the insulation layer (from "Launch and Entry Vehicle Design, Univ. Maryland, D.L. Akin")
M.tank_lox = rho * S_lox * t_lox + M.lox_insulation; %mass of the empty lox tank [kg]
M.tot_lox = M.tank_lox + mlox; %[kg] mass of full lox tank

%last connector (between second tank and motors)
shape3.r = R_lox;
shape3.h = (2/3) * R_lox + h_dome_lox; 
load3 = loads;
load3.m = m_sust + C1.m + M.tot_rp1 + C2.m + M.tot_lox; %accounts for sustained masses : upper stages, first and second connector, fairing, rp1 tank and lox tank structural mass (approximated as 11% of lox mass)
load3.F_drag = loads.F_drag; %aerodynamic force [N]
if nargin < 8
    [C3.m, C3.th] = buckling_bending(shape3, load3, mat, 0);
else
    shape3.h0 = h0_C3;
    [C3.m, C3.th, C3.XY] = buckling_bending(shape3, load3, mat, 0, 1);
    plot(C3.XY(1,:), C3.XY(2,:), '--k', DisplayName='true');
end
h.CG.C3 = h0_C3 + shape3.h / 2;%[m]

%TOTAL MASSES:
M.tot = M.tot_lox + M.tot_rp1 + M.motor + M.fairing + M.avionics + M.T_struct + M.wiring + C1.m + C2.m + C3.m; %[kg] total mass of motors, tanks, propellant, connectors, wiring, avionics, structures
M.tanks = M.tank_lox + M.tank_rp1; %[kg] mass of the two empty tanks
M.str = M.tanks + M.motor + M.fairing + M.avionics + M.T_struct + M.wiring + C1.m + C2.m + C3.m; %[kg] inert mass of stage (motors, tanks, connectors, stage fairing, wiring, avionics, structures)
M.eps = M.str / M.tot;
M.C1 = C1;
M.C2 = C2;
M.C3 = C3;

%HEIGHTS
h.tank_lox = l_lox(t_lox); %[m] height of tank
h.tank_rp1 = l_rp1(t_rp1); %[m] height of tank
h.C1 = shape1.h; %[m] first connector / top interstage
h.C2 = shape2.h; %[m] second connector / intertank
h.C3 = shape3.h; %[m] thirk connector / aft skirt
h.tot = h_cyl_lox + h_cyl_rp1 + h.C2 + h.C3; %[m] total height of stage
h.CG.avionics = h0_C1 + R_rp1; %[m]
h.CG.T_struct = h_motor + 0.5 * R_lox; %[m]
h.CG.tot = (h.CG.lox * M.tot_lox + h.CG.rp1 * M.tot_rp1 +...
        + h.CG.C1 * C1.m + h.CG.C2 * C2.m + h.CG.C3 * C3.m +...
        + h.CG.avionics * M.avionics + 0.5 * h_motor * M.motor +...
        + h.CG.T_struct * M.T_struct) / M.tot;

%PLOT OF STAGE:
if nargin > 7
    if R_lox == R_rp1
        c = 1;
    else
        c = -1;
    end
        %lox
        bottom = @(k) h0_lox + h_dome_lox -sqrt(R_lox^2 - k.^2)/AR_lox;
        top = @(k) bottom(R_lox) + sqrt(R_lox^2 - k.^2)/AR_lox + h_cyl_lox;
        K = linspace(-R_lox, R_lox, 1e4);
        plot(K, bottom(K), 'k'); grid on; axis equal; hold on;
        plot(K, top(K), 'k');
        plot([ R_lox,  R_lox], [bottom(R_lox), bottom(R_lox) + h_cyl_lox], 'k');
        plot([-R_lox, -R_lox], [bottom(R_lox), bottom(R_lox) + h_cyl_lox], 'k');
        
        %rp1 
        bottom = @(k) h0_rp1 + c* sqrt(R_rp1^2 - k.^2)/AR_rp1;
        top = @(k) bottom(R_rp1) + sqrt(R_rp1^2 - k.^2)/AR_rp1 + h_cyl_rp1;
        K = linspace(-R_rp1, R_rp1, 1e4);
        plot(K, bottom(K), 'k');
        plot(K, top(K), 'k');
        plot([ R_rp1,  R_rp1], [bottom(R_rp1), bottom(R_rp1) + h_cyl_rp1], 'k');
        plot([-R_rp1, -R_rp1], [bottom(R_rp1), bottom(R_rp1) + h_cyl_rp1], 'k');

    if h_cyl_lox == 0 %sphere
        plot(0, h0_lox + h_dome_lox, '+k');
    end
    if h_cyl_rp1 == 0 %sphere
        plot(0, h0_rp1 + h_dome_rp1, '+k');
    end
    xlabel('x [m]', 'Interpreter','latex');
    ylabel('y [m]', 'Interpreter','latex');
end
end

function [m, th, XY] = buckling_bending(shape, load, mat, plotcase)

% based on NASA papers in shared folder (SP-8007 & SP-8019)
% computes connectors masses, heights and thicknesses to sustain
% compression loads and avoid buckling effect

%constants:
g = 9.81; %[m/s^2] gravitational acceleration

%recover loads:
m = load.m; %sustained mass [kg]
n = load.n; %longitudinal load factor [-]
K = load.K; %factor of safety [-]
F_aero = load.F_drag; %aerodynamic drag force [N]
p = load.p; %internal pressure [Pa]
M_exp = load.M_exp; %expected flessional 

%recover material characteristics:
id = mat.ID; %[-] ID of the material: 1 for Ti, 2 for Al 2XXX, 3 for Steel, 4 for Carbon Fiber, 5 for Al 7XXX
E = mat.E; %[Pa] Young modulus
rho = mat.rho; %[kg/m^3]
t_min = mat.t_min; %[m] minimum thickness for manufacturability 
sy = mat.sy; %[Pa] tensile yield stress
su = mat.su; %[Pa] tensile ultimate stress
nu = mat.nu; %[-] Poisson's ratio

%recover dimensions:
r = shape.r;
h = shape.h; %distance between the base of the domes
if length(r) > 1 %trucated-cone
    %recover cone shape characteristics
    alpha = asin( ( r(2) - r(1) ) / h );
    L = sqrt( h^2 - (r(2) - r(1))^2 );
    l = cos(alpha) * L; %height of the shell
    rho1 = r(1);
    rho2 = r(2);
    r(2) = cos(alpha) * r(2); 
    r(1) = cos(alpha) * r(1);

    %plot
    if nargin > 4
        h0 = shape.h0; %height at which the connector is placed 
        %for the plotting
        y = h0 + [rho2*sin(alpha/2), rho2*sin(alpha/2), rho2*sin(alpha/2)+l, rho2*sin(alpha/2)+l, rho2*sin(alpha/2), rho2*sin(alpha/2)];
        XY = [0, rho2, rho1, -rho1, -rho2, 0; y];
    end
else
    %recover cylinder shape characteristics
    alpha = 0;
    L = h;
    if nargin > 4
        h0 = shape.h0; %height at which the connector is placed 
        %for the plotting
        y = h0 + [0, 0, L, L, 0, 0];
        XY = [0, r, r, -r, -r, 0; y];
    end
end

%compute sustained load:
F_load = m * n * g + F_aero; %load [N]

%get the wall flexural stiffness per unit width:
D = @(th) E * th^3 / ( 12 * (1-nu^2) );

%compute delta_gamma for pressure-increased performances:
dg = @(th) d_gamma(p, E, r, th, alpha); 

%6 SITUATIONS: p OR NON p , CONICAL/CYLINDRICAL (BUT THERE AREN'T CONICAL TANKS), ISOTROPIC/ORTHOTROPIC
switch id
    case 4 %(for CF, orthotropic expressions)

    otherwise
        if alpha == 0 %(cylindrical shape, NASA SP-8007)
            k1 = 0.8;

            %compute knockdown factor:
            phi = @(th) (1/16) * sqrt(r(1)/th);
            gP = @(th) 1 - 0.901 * ( 1 - exp( -phi(th) ) );
            gM = @(th) 1 - 0.731 * ( 1 - exp( -phi(th) ) );

            %in cyl, distinguish between presurrized and unpressurized cases:
            if p == 0
                k2 = 1;
            else %( p > 0 )
                k2 = 0;
            end
        else %(alpha > 0) (conical shape, NASA SP-8019)
            k1 = 0.5;
            k2 = 0;

            %compute knockdown factor:
            gP = @(th) 0.33; 
            gM = @(th) 0.41;
        end
end

%compute critical loads expressions
Pcr = @(th) k2 * k_x(nu, L, r, th, gP) * 2 * pi^3 * D(th) * r   / L^2 +...   %if cyl and p=0 (only metals)
         + (1-k2) * ( 2*pi    * E * th^2 * ( gP(th) / sqrt( 3 * (1-nu^2) ) + dg(th) ) + p * pi * r(1)^2 ); %in any other case (only metals)
Mcr = @(th) k2 * k_x(nu, L, r, th, gM)   *   pi^3 * D(th) * r^2 / L^2 +...   %if cyl and p=0 (only metals)
         + (1-k2) * pi*r(1) * E * th^2 * ( gM(th) / sqrt( 3 * (1-nu^2) ) + dg(th) ) + p * pi * r(1)^2 * k1;%in any other case (only metals)

%relation to be satisfied in combined stress condition
f = @(th) K*F_load/Pcr(th) + K*M_exp/Mcr(th) - 1; %this must be <0 to preserve the structure

%minimum thickness to sustain only Pcr:
th = fzero( f , [0, 1]); %[m]

%check on manufacturability (th>t_min) and on validity of the NASA model (th>r/1500):
th = max(th, t_min, r/1500); %[m] 

%compute the mass
if length(r) > 1
    S = pi * L * ( r(2) + r(1) ); %surface of the truncated cone
    m = S * th * rho; %[kg]
else
    m = 2*pi * r * h * th * rho; %[kg]
end
end

function dg = d_gamma(p, E, r, th, alpha)

% based on NASA paper SP-8007-2020/REV 2:

if p == 0
    dg = 0;
else
    if nargin < 5
        alpha = 0;
    end
    param = (p/E)*(r(1)/ (th*cos(alpha)) )^2;

    % from curve fitting with 
    % param = [0.02, 0.1, 1, 10]
    % d_gamma = [0.026, 0.09, 0.2, 0.23]
    % and f(x) = (a*x+b)/(c*x+d)
    % we obtained the fitting (R-square = 0.99995, SSE = 1.241e-6, DFE = 0):

    a  = 0.3045;   
    b  = 0.0001;
    c  = 1.3064;
    d  = 0.2104;

    dg = ( a * param + b ) / ( c * param + d );
end
end

function kx = k_x(nu, L, r, th, g)

% based on NASA paper SP-8007-2020/REV 2:

gZ = @(th) g * L^2 / (r*th) * sqrt(1-nu^2);

% from curve fitting with 
% param = [0.1, 1, 3, 100, 1000]
% d_gamma = [1, 1.1, 2, 64, 700]
% and f(x) = c*(a+x^2)^b
% we obtained the fitting (R-square = 1.0000, SSE = 0.0025, DFE = 1):
 
a = 3.2513;    
b = 0.5195;   
c = 0.5348;

kx = c * (a + gZ(th)^2 ) ^ b;
end

function [fairing] = fairing_fun(m_pay_max, m_pay, fairing, loads, plotcase)

% computes the shape and mass of the fairing.
% based on cubesats volumetric density constraints, takes as input also the
% aerodynamic loads.
% fairing is simplified as a cylinder and a cone.

%cubesats characteristics:
vol_den = 1.33 * 1e3; %[kg/m^3] volumetric density of cubesats
a = 0.1; %[m] cubesat unit edge length 
b = 0.05; %[m] margin distance between payload and fairing

%recover fixed shape characteristics
d0 = fairing.base_diam; %[m] diameter of the base
L_nose = fairing.nose_length; %[m] length of the conical nose

%recover loads:
n = loads.n; %longitudinal load factor [-]
K = loads.K; %factor of safety [-]

%compute payload volume (considering cubesats sizes)
V = m_pay_max / vol_den; %[m^3] volume of max payload
V_real = m_pay / vol_den; %[m^3] volume of real payload

%find usable base for payload
y = 0 : a : d0/2; 
N = length(y) - 1; 
base = 0;
x = zeros(N, 1);
for j = 1 : N
    x(j) = floor( sqrt( 0.25*d0^2 - y(j)^2 - b ) / a ) * a; %[m] at y(i) level, how much space i have to fill with cubesats
    base = base + 4 * a * x(j); %[m^2] update base value 
end
L_min = V / base; %[m] height of the ammissible payload
L_real = V_real / base; %[m] height of the real payload 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get adapter characteristics
adapter_loads.m = m_pay; %sustained mass [kg]
adapter_loads.n = n; %longitudinal load factor [-]
adapter_loads.K = K; %factor of safety [-]
adapter_loads.F_drag = 0; %aerodynamic drag force [N] is null for the adapter
adapter.mat_id = fairing.mat_id; %in reality, we just need the height of the adapter for this step, therefore we can put any material, since the shape is fixed
adapter.base_diam = fairing.base_diam; %[m] diameter of the base
adapter = adapter_fun(adapter, adapter_loads);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find height of the cylindrical part of the fairing
L1 = L_min + b + adapter.h; %[m] 

%surface of the fairing
S_cyl = L1 * pi * d0; %[m^2] cylinder surface
S_cone = pi * sqrt( d0^2 / 4 + L_nose^2 ) * d0/2; %[m^2] conical surface
fairing.S = S_cyl + S_cone; %[m^2] total surface

%compute mass (from "Launch and Entry Vehicle Design, Univ. Maryland, D.L.Akin")
fairing.m = 4.95 * fairing.S ^ 1.15; %[kg]  

%retrieve fairing parameters
% fairing.M = nose.M + cyl.M; %[kg] total mass
fairing.L = L1 + L_nose; %[m] total height
fairing.d = d0; %[m] maximum diameter
% fairing.th_cone = nose.th; %[m] thickness of the nose
% fairing.th_cyl = cyl.th; %[m] thickness of the cylinder
fairing.V = V; %[m^3] maximum volume for the payload

%plot of the fairing
if nargin > 4
    h0 = fairing.h0;
    X = [d0/2, d0/2, 0, -d0/2, -d0/2];
    Y = h0 + [0, L1, fairing.L, L1, 0];
    plot(X, Y, 'k', DisplayName='true'); grid on; hold on; axis equal;%fairing
    X_a = [0, adapter.r(2), adapter.r(1), -adapter.r(1), -adapter.r(2), 0];
    Y_a = h0 + [0, 0, adapter.h, adapter.h, 0, 0];
    plot(X_a, Y_a, '--r', DisplayName='true'); %adapter
    X_p = [0, x(1), x(1), -x(1), -x(1), 0];
    Y_p = h0 + [adapter.h, adapter.h, adapter.h + L_real, adapter.h + L_real, adapter.h, adapter.h]; 
    plot(X_p, Y_p, 'b', DisplayName='true'); %payload
    X_r = [0, 0.6, 0.6, -0.6, -0.6, 0];
    Y_r = [0, 0, h0, h0, 0, 0];
    plot(X_r, Y_r, '--k', DisplayName='true');
    legend('Fairing', 'Adapter', 'Payload','Hypothetical Rocket', 'interpreter', 'latex');
    xlabel('x [m]', 'Interpreter','latex');
    ylabel('y [m]', 'Interpreter','latex');
end

end

function [adapter] = adapter_fun(adapter, loads)

% this function computer adapter characteristics.
% it uses the loads, the fixed maximum diameter, adn the material

%recover fixed shape characteristics
d0 = adapter.base_diam; %[m] diameter of the base
adapter.r = [d0/3, d0/2]; %[m] payload adapter r1, r2
adapter.h = d0 / 4; %[m] payload adapter height 

%update internal pressure load:
loads.p = 0; %[Pa]

%recover material properties
mat = mat_switch(adapter.mat_id);

%compute adapter characteristics
[adapter.m , adapter.th] = buckling_bending(adapter, loads, mat);

end

function [mat] = mat_switch(mat_id)

%this function gets the info about the selected material:
% mat_id = 1 for Ti, 2 for Al 2XXX, 3 for Steel, 4 for Carbon Fiber, 5 for Al 7XXX

%switch
switch mat_id 
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

%recover material properties:
mat.ID = mat_id; %[-] ID of the material: 1 for Ti, 2 for Al 2XXX, 3 for Steel, 4 for Carbon Fiber, 5 for Al 7XXX
mat.rho = rho; %[kg/m^3] material density
mat.t_min = t_min; %[m] material manufacturability minimum thickness
mat.E = E; %[Pa] Young modulus
mat.sy = sy; %[Pa] yelding stress
mat.su = su; %[Pa] ultimate stress
mat.nu = nu; %[-] Poisson's ratio

end


