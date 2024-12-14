

%% (only for flight) run this with the matlab structs with the data 
clc;
close all;
load('T.mat');
load('parout.mat');
load('Y.mat');

%row:
r = 1;

%earth radius:
Re = 6378000; %[m]

%fos:
FoS = loads_it(r).K;
%FoS = 0.5;

%materials:
mat1 = M_it(r).mat1;
mat2 = M_it(r).mat2;

%get data:
h_tot = h_it(r).tot;
h_com = h_it(r).CG;
h_st1 = h_it(r).stg1.tot;
h_tot = h_tot - h_it(r).stg1.motor - (2/3)*h_it(r).fairing;
h_com = h_com - h_it(r).stg1.motor;
h_st1 = h_st1 - h_it(r).stg1.motor;

%diameters:
d2 = M_it(r).diam2;
d1 = M_it(r).diam1;

%cross area:
S2 = pi * d2^2 / 4;
S1 = pi * d1^2 / 4;

%get thicknesses of second stage:
t2(1) = M_it(r).th2.C1;
t2(2) = M_it(r).th2.rp1;
t2(3) = M_it(r).th2.lox;
t2(4) = M_it(r).th2.C3;
%get thicknesses of first stage:
t1(1) = M_it(r).th1.C1;
t1(2) = M_it(r).th1.rp1;
t1(3) = M_it(r).th1.lox;
t1(4) = M_it(r).th1.C3;
%get lenghts of components of second stage:
L2(1) = h_it(r).stg2.C1;
L2(2) = h_it(r).stg2.cyl_rp1;
L2(3) = h_it(r).stg2.cyl_lox;
L2(4) = h_it(r).stg2.C3;
%get lenghts of components of second stage:
L1(1) = h_it(r).stg1.C1;
L1(2) = h_it(r).stg1.cyl_rp1;
L1(3) = h_it(r).stg1.cyl_lox;
L1(4) = h_it(r).stg1.C3;
%internal pressures of second stage:
p2(1:4) = 500*1e3;
p2(1) = 0;
p2(4) = 0;
%internal pressures of first stage:
p1 = p2;
%get sustained masses of second stage:
m2(1) = M_it(r).fairing + M_it(r).pay_effective + M_it(r).stg2.avionics + M_it(r).stg2.wiring;
m2(2) = m2(1) + M_it(r).stg2.C1.m;
m2(3) = m2(2) + M_it(r).stg2.tot_rp1;
m2(4) = m2(3) + M_it(r).stg2.tot_lox;
%get sustained masses of second stage:
m1(1) = M_it(r).stg2.tot + M_it(r).stg1.avionics + M_it(r).stg1.wiring + m2(1);
m1(2) = m1(1) + M_it(r).stg1.C1.m;
m1(3) = m1(2) + M_it(r).stg1.tot_rp1;
m1(4) = m1(3) + M_it(r).stg1.tot_lox;

%correlate with nx, ny:
R = sqrt(Y(:,1).^2 + Y(:,2).^2);
ng = -(Re./R).^2;
alp = atan(Y(:,1)./Y(:,2));
bet = atan(-Y(:,4)./Y(:,3));
gamma = alp - bet; %[rad] flight path angle
vt = cos(gamma) .* sqrt( Y(:,4).^2 + Y(:,3).^2 ); %[m/s] transversal velocity (in orbit ref frame)
ncp = vt.^2 ./ (9.81*R);
ndiff = ng + ncp; %[-] local vertical load factor
eps = gamma + parout.alpha; %[rad] LVLH pitch 

nx_old = parout.acc(:, 1)/9.81; %[-] longitudinal old
nz_old = parout.acc(:, 2)/9.81; %[-] transversal old
nx_t = nx_old - ncp .*cos(eps); %[-] long with centrifugal force
nz_t = nz_old + ncp .*sin(eps); %[-] tran with centrifugal force

Qdyn = parout.qdyn; %[Pa] dynamic pressure
M_t = parout.m_vec; %[kg] mass of the rocket 
m = length(nx_t);

%P load for second stage:
% m = 1e3;
% P2_min = 0;
% P2_max = loads_it(r).F_drag(3) + M_it(r).nx2 * M_it(r).M1;
% P2 = linspace(P2_min, P2_max, m);
alpha = parout.alpha;
cd2 = parout.coeffs(:,1);
cl2 = parout.coeffs(:,2);
cx2 = cd2 .* cos(alpha) - cl2 .* sin(alpha); 
P2 = zeros(m, 4);
for j = 1:4
    P2(:, j) = Qdyn' * S2 .* cx2' + nx_t' * 9.81 .* m2(j);
end

%P load for first stage:
% P1_min = 0;
% P1_max = loads_it(r).F_drag_tot + M_it(r).nx1 * M_it(r).M0;
% P1 = linspace(P1_min, P1_max, m);
cd1 = parout.coeffs(:,1);
cl1 = parout.coeffs(:,2);
cx1 = cd1 .* cos(alpha) - cl1 .* sin(alpha); 
P1 = zeros(m, 4);
for j = 1:4
    P1(:, j) = Qdyn' * S1 .* cx1' + nx_t' * 9.81 .* m1(j);
end

%actual bending moment due to nz:
Mz1 = abs(nz_t .* M_t * 9.81 * h_com * ( 1 - h_com/h_tot ) );
Mz2 = Mz1 * ( h_tot - h_st1 ) / ( h_tot - h_com );


%SECOND stage maximum bending moment:
M2maxj = zeros(m, 4);%initialize
for j = 1:4 %each column of of M2maxj represents the maximum bending moment that the j-th component can take for the given axial compressive load P
    M2maxj(:, j) = admissible_M_in_flight(diam2, t2(j), p2(j), FoS, mat2, P2(:, j), L2(j));
end
M2max = zeros(m, 1);%initialize
for j = 1:m %select only the minimum admissible value
    M2max(j) = min(M2maxj(j, :));
end

%FIRST stage maximum bending moment:
M1maxj = zeros(m, 4);%initialize
for j = 1:4 %each column of of M2maxj represents the maximum bending moment that the j-th component can take for the given axial compressive load P
    M1maxj(:, j) = admissible_M_in_flight(diam1, t1(j), p1(j), FoS, mat1, P1(:, j), L1(j));
end
M1max = zeros(m, 1);%initialize
for j = 1:m %select only the minimum admissible value
    M1max(j) = min(M1maxj(j, :));
end

% figure(1)
% plot( P1, M1max, '--r'); grid on; axis equal; hold on;
% xlabel('Compressive Axial Load P [N]', 'Interpreter','latex');
% ylabel('Admissible Bending Moment M [Nm]', 'Interpreter','latex');
% title('Admissible for $1^{st}$ Stage', 'Interpreter','latex');
% 
% figure(2)
% plot( P2, M2max, '--b'); grid on; axis equal;
% xlabel('Compressive Axial Load P [N]', 'Interpreter','latex');
% ylabel('Admissible Bending Moment M [Nm]', 'Interpreter','latex');
% title('Admissible for $2^{nd}$ Stage', 'Interpreter','latex');

figure(1)
plot( T(1:3610), M1max(1:3610), '-r'); grid on; hold on;%axis equal; hold on;
plot( T(1:3610), Mz1(1:3610), '-k'); 
xlabel('Time in flight [s]', 'Interpreter','latex');
ylabel('Bending Moment M [Nm]', 'Interpreter','latex');
title('Admissible for $1^{st}$ Stage', 'Interpreter','latex');
legend('Admissible', 'Effective', 'Interpreter','latex');

figure(2)
plot( T, M2max, '-b'); grid on; hold on;%axis equal;
plot( T, Mz2, '-k');
xlabel('Time in flight [s]', 'Interpreter','latex');
ylabel('Bending Moment M [Nm]', 'Interpreter','latex');
title('Admissible for $2^{nd}$ Stage', 'Interpreter','latex');
legend('Admissible', 'Effective', 'Interpreter','latex');


%% Functions

function [Mmax] = admissible_M_in_flight(d, th, p, FoS, mat, P, L)
% for a given axial compressive load P, computes the maximum bending Moment
% Mmax

% ONLY FO CYLINDRICAL VESSELS, AND FOR p_hydro ==0:
alpha = 0;
k1 = 0.8;
p_hydro = 0;

%get number of iterations:
n = length(P);

%radius of rocket:
r = d/2;

%material properties:
mat = mat_switch(mat);
E = mat.E;
nu = mat.nu;
sy = mat.sy;

%phi:
phi = (1/16)*sqrt(r/th);

%gp, gm, dg:
gP = 1 - 0.901 * ( 1 - exp( -phi ) );
gM = 1 - 0.731 * ( 1 - exp( -phi ) );
dg = d_gamma(p, E, r, th, alpha);

%kx:
kx = k_x(nu, L, min(r), th, gP);

%D:
D = E * th.^3 / ( 12 * (1-nu^2) );

%initialize:
Mmax = zeros(n, 1);

for i = 1:n

    P_actual = P(i);

    %buckling:
    if p == 0
        Pcr = kx * 2 * pi^3 * D * min(r)   / L^2;
        Mcr = kx   *   pi^3 * D * min(r)^2 / L^2;
    else
        Pcr = ( 2*pi    * E * th.^2 .* ( gP ./ sqrt( 3 * (1-nu^2) ) + dg ) + p * pi * min(r)^2 );
        Mcr = pi*min(r) * E * th.^2 .* ( gM ./ sqrt( 3 * (1-nu^2) ) + dg ) + p * pi * min(r)^2 * k1;
    end
    Mb = Mcr * (1/FoS - P_actual/Pcr);

    %yielding:
    My = pi * r^2 * (th*sy/FoS - P_actual/(2*pi*r) + (p+p_hydro)*r/2);
    
    %admissible for this P load
    Mmax(i) = min(Mb, My);
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

gZ = @(t) g * L^2 / (r*t) * sqrt(1-nu^2);

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
    case 4 % Carbon fiber  Toray M46J
        rho = 1600; %[kg/m^3]
        t_min = 0.90 * 1e-3; %[m] minimum thickness for manufacturability
        E = 222 * 1e9; %[Pa] young modulus
        sy = 1090 * 1e6; %[Pa] tensile yield stress
        su = sy; %[Pa] tensile ultimate stress
        nu = 0.28; %[-] Poisson's ratio
    case 5 % Al 7075 T6
        rho = 2810; %[kg/m^3]
        t_min = 1.06 * 1e-3; %[m] minimum thickness for manufacturability
        E = 70 * 1e9; %[Pa] young modulus
        sy = 440 * 1e6; %[Pa] tensile yield stress
        su = 517 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.33; %[-] Poisson's ratio
    case 6 % AlLi (2090)
        rho = 2590; %[kg/m^3]
        t_min = 0.5 * 1e-3; %[m] minimum thickness for manufacturability  
        E = 76 * 1e9; %[Pa] young modulus
        sy = 500 * 1e6; %[Pa] tensile yield stress
        su = 550 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.34; %[-] Poisson's ratio
    case 7 % Carbon fiber Hexcel® HexTow® IM7
        rho = 1600; %[kg/m^3]
        t_min = 0.90 * 1e-3; %[m] minimum thickness for manufacturability
        E = 148 * 1e9; %[Pa] young modulus
        sy = 1000 * 1e6; %[Pa] tensile yield stress
        su = sy; %[Pa] tensile ultimate stress
        nu = 0.28; %[-] Poisson's ratio
    case 8 % Al 6061 T6
        rho = 2700; %[kg/m^3]
        t_min = 0.5 * 1e-3; %[m] minimum thickness for manufacturability
        E = 68.9 * 1e9; %[Pa] young modulus
        sy = 276 * 1e6; %[Pa] tensile yield stress
        su = 310 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.33; %[-] Poisson's ratio
    case 9 % 300M Steel alloy
        rho = 7830; %[kg/m^3]
        t_min = 0.25 * 1e-3; %[m] minimum thickness for manufacturability
        E = 207 * 1e9; %[Pa] young modulus
        sy = 1586 * 1e6; %[Pa] tensile yield stress
        su = 1931 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.3; %[-] Poisson's ratio
    case 10 % 2219 Al alloy
        rho = 2840; %[kg/m^3]
        t_min = 0.5 * 1e-3; %[m] minimum thickness for manufacturability
        E = 73.1 * 1e9; %[Pa] young modulus
        sy = 350 * 1e6; %[Pa] tensile yield stress
        su = 440 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.33; %[-] Poisson's ratio
    case 11 % Carbon fiber  (MatWeb)
        rho = 1420; %[kg/m^3]
        t_min = 3 * 1e-3; %[m] minimum thickness for manufacturability
        E = 101 * 1e9; %[Pa] young modulus
        sy = 1260 * 1e6; %[Pa] tensile yield stress
        su = sy; %[Pa] tensile ultimate stress
        nu = 0.286; %[-] Poisson's ratio
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
