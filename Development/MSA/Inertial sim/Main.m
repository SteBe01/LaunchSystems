%% Test

clear, clc
close all

param.T = 254e3;
param.d = 1.8;
param.S = param.d^2/4*pi;
param.Cd = 0.7;
param.T_direction = deg2rad(30);
param.m_MECO = 1000;

param.g0 = 9.81;
param.Isp = 280;
param.m_dot = -param.T/(param.Isp*param.g0);

% initial conditions
x0 = 0;
y0 = 12000;
vx0 = 0;
vy0 = 350;
m0 = 11500;
y0 = [x0 y0 vx0 vy0 m0]';

tf = 1e5;

odefun = @(t,x) dynamics(t,x,param);
options = odeset('RelTol',1e-8, 'AbsTol',1e-8, 'Events', @(t,y) touchdown(t,y));
[t,y] = ode89(odefun,[0 tf],y0,options);

figure, hold on, grid on, axis equal, xlabel("x [km]"), ylabel("y [km]")
plot(y(:,1).*1e-3,y(:,2).*1e-3)

figure
subplot(2,2,1)
plot(t,y(:,1).*1e-3), grid on, xlabel("Time [s]"), ylabel("x [km]")
subplot(2,2,2)
plot(t,y(:,2).*1e-3), grid on, xlabel("Time [s]"), ylabel("y [km]")
subplot(2,2,3)
plot(t,y(:,3).*1e-3), grid on, xlabel("Time [s]"), ylabel("vx [km/s]")
subplot(2,2,4)
plot(t,y(:,4).*1e-3), grid on, xlabel("Time [s]"), ylabel("vy [km/s]")

figure, hold on, grid on, xlabel("Time [s]"), ylabel("Mass [kg]")
plot(t,y(:,5))

H0 = 4;
rho0 = 1.225;
rho = @(h) rho0*exp(-h.*1e-3/H0);
dynPress = @(h,v) 0.5*rho(h)*v^2;
test = zeros(length(y),1);
for ii = 1:length(y)
    test(ii) = dynPress(norm(y(ii,1),y(ii,2)),norm(y(ii,3),y(ii,4)));
end

figure, hold on, grid on, xlabel("Time [s]"), ylabel("Dynamic pressure [kPa]")
plot(t,test.*1e-3)

test = zeros(length(y),1);
for ii = 1:length(y)
    test(ii) = rad2deg(atan2(y(ii,4),y(ii,3)));
end

figure, hold on, grid on, xlabel("Time [s]"), ylabel("Pitch angle [deg]")
plot(t,test)



%% Functions

function [value, isterminal, direction] = touchdown(~, y)
    value = y(2);
    isterminal = 1;
    direction = -1;
end

function dx = dynamics(~,x,param)

% x(1) = x;
% x(2) = y;
% x(3) = vx;
% x(4) = vy;
% x(5) = m;

T = param.T;
S = param.S;
Cd = param.Cd;

% air density
h0_vect = [0 25 30 40 50 60 70 80 90 100 110 120 130 140 150 ...
    180 200 250 300 350 400 450 500 600 700 800 900 1000]';
rho0_vect = [1.225 3.899*1e-2 1.774*1e-2 3.972*1e-3 1.057*1e-3 ...
    3.206*1e-4 8.770*1e-5 1.905*1e-5 3.396*1e-6 5.297*1e-7 9.661*1e-8 ...
    2.438*1e-8 8.484*1e-9 3.845*1e-9 2.070*1e-9 5.464*1e-10 ...
    2.789*1e-10 7.248*1e-11 2.418*1e-11 9.158*1e-12 3.725*1e-12 ...
    1.585*1e-12 6.967*1e-13 1.454*1e-13 3.614*1e-14 1.170*1e-14 ...
    5.245*1e-15 3.019*1e-15]';
H_vect = [7.249 6.349 6.682 7.554 8.382 7.714 6.549 5.799 5.382 5.877 ...
    7.263 9.473 12.636 16.149 22.523 29.740 37.105 45.546 53.628 ...
    53.298 58.515 60.828 63.822 71.835 88.667 124.64 181.05 268.00]';

height = x(2).*1e3;
idx = sum(h0_vect < height) + 1;
if idx < length(h0_vect)
    rho0 = rho0_vect(idx);
    h0 = h0_vect(idx);                            % [km]
    H = H_vect(idx);                              % [km]
    
    rho = rho0 * exp(-(height-h0)/H);                % [kg/m^3]
else
    rho = 0;
end

D = @(v) 0.5*rho*v^2*S*Cd;

T_direction = param.T_direction;
Tx = cos(T_direction)*T;
Ty = sin(T_direction)*T;

Dx = D(x(3));
Dy = D(x(4));

m = x(5);

MECO = 0;
if m < param.m_MECO
    MECO = 1;
end

if ~MECO
    dx = zeros(5,1);
    dx(1:2) = x(3:4);
    dx(3) = Tx/m - Dx/m;
    dx(4) = Ty/m - Dy/m - param.g0;
    dx(5) = param.m_dot;
else
    dx = zeros(5,1);
    dx(1:2) = x(3:4);
    dx(3) = - Dx/m * sign(x(3));
    dx(4) = - Dy/m * sign(x(4)) - param.g0;
    dx(5) = 0;
end

end

