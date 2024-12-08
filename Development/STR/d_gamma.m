
clc;
clear;
close all;

%material
mat_id = 1; %Ti

%material properties
mat = mat_switch(mat_id);

%thickness
th = 3e-3; %[m]

%radius
r = 0.6; %[m]

%parameter
param = @(p) p./mat.E * (r/th)^2;%[-]

%how the param varies with pressure
p = [0.3, 2.8, 50]*1e6;
figure(1);
loglog(p, param(p), '+');grid on;

%delta_gamma
%hypothesis
h = 1.5*1e1;
dg = @(p) 0.2 * h * param(p) ./ ( h * param(p) + 3 ) + 0.2 / ( 10 * h + 1 );
p = linspace(0.003, 50, 1e3) * 1e6;
figure(2);
loglog(param(p), dg(p));grid on;hold on;
loglog([0.02, 0.1, 1, 10], [0.02, 0.09, 0.2, 0.2]);

%f = fit([0.02, 0.1, 1, 10], [0.02, 0.09, 0.2, 0.2], 'rat23');
% curveFitter([0.02, 0.1, 1, 10], [0.026, 0.09, 0.2, 0.23]);




%%
clc;
clear;
close all;

%material
mat_id = 1; %Ti

%material properties
mat = mat_switch(mat_id);

%thickness
th = 3e-3; %[m]

%radius
r = 0.6; %[m]

%parameter
param = @(p) p./mat.E * (r/th)^2;%[-]

% from curve fitting with 
% param = [0.02, 0.1, 1, 10]
% d_gamma = [0.026, 0.09, 0.2, 0.23]
% and f(x) = (a*x+b)/(c*x+d)
% we obtained the fitting (R-square = 0.99995, SSE = 1.241e-6, DFE = 0):

a  = 0.3045;   
b  = 0.0001;
c  = 1.3064;
d  = 0.2104;

dg = @(p) ( a * param(p) + b )./( c * param(p) + d );

%pressures
p = linspace(0.003, 50, 1e3) * 1e6;
figure(3);
loglog(param(p), dg(p));grid on;hold on;

%new plot
d_g = @(param) ( a * param + b )./( c * param + d );
param = linspace(0.02, 100, 1e6);
figure(4);
loglog(param, d_g(param), '.');grid on;axis equal;hold on;
ylabel('$\Delta\gamma$', 'Interpreter','latex');
xlabel('$\frac{p}{E}\left (\frac{r_1}{t \cos{\alpha}}\right )^2$','Interpreter','latex');


%retake to confirm that our parameters fall in the casistics
param = @(p, th) p./mat.E * (r/th)^2;%[-]
d_g = @(p, th) ( a * param(p, th) + b )./( c * param(p, th) + d );
p = [0.3, 2.8, 50]*1e6;
th = [0.25, 1, 1.5, 3]*1e-3;

n = length(p);
m = length(th);

for i = 1 : n
    for j = 1 : m
        prm = param(p(i), th(j));
        dgm = d_g(p(i), th(j));
        loglog(prm, dgm, '+'); hold on;
    end
end

%%

figure(5)
surf(param(p, th), p, th);

%%
x = linspace(0, 1000, 1e4);
y = linspace(0, 5000, 1e4);
figure(6)
loglog(x, y);

%%
%curveFitter([0.1, 1, 3, 100, 1000], [1, 1.1, 2, 64, 700]);

% from curve fitting with 
% param = [0.1, 1, 3, 100, 1000]
% d_gamma = [1, 1.1, 2, 64, 700]
% and f(x) = c*(a+x^2)^b
% we obtained the fitting (R-square = 1.0000, SSE = 0.0025, DFE = 1):
 
a = 3.2513;    
b = 0.5195;   
c = 0.5348;

gZ = linspace(0.1, 1e4, 1e5);

f = @(x) c*(a+x.^2).^b;

loglog(gZ, f(gZ));grid on;

%%




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

