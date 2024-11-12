%% Cannon shooting

clear, clc
close all

input.vmag = 30;
input.theta = deg2rad(45);
input.tflight = 50;

target.x = 1000;
target.y = 1000;

param.g = 9.81; 
param.nGrid = 300; 


%% Optimization

traj = traj_fun(input, param);

guess.dx0 = traj.dx(1);  %Initial horizontal speed
guess.dy0 = traj.dy(1);  %Initial vertical speed
guess.t   = traj.t(end); %Time of flight
problem.x0    = [guess.dx0; guess.dy0; guess.t];

problem.Aineq = [];
problem.Bineq = [];
problem.Aeq   = [];
problem.Beq   = [];
problem.lb    = [];
problem.ub    = [];

problem.solver = 'fmincon';
problem.options = optimset('Display','iter','MaxIter',2000);

problem.objective = @(var_input) obj_fun(var_input(1), var_input(2));  %Objective (cost) function
problem.nonlcon  = @(var_input) nl_cons(var_input, target, param);     %NonLinear constraints

[xSoln, fVal, exitFlag] = fmincon(problem);

final_sol.vmag = sqrt(xSoln(1)^2+xSoln(2)^2);
final_sol.theta = atan2(xSoln(2),xSoln(1));
final_sol.tflight = xSoln(3);
final_traj = traj_fun(final_sol, param);

figure, hold on, grid on
plot(final_traj.x,final_traj.y)






%% Functions

function cost = obj_fun(dx, dy)

m = 1;                              % [kg] - to be changed for different masses
cost = 0.5*m*(dx.*dx+dy.*dy);

end

function sol = traj_fun(input, param)

v0 = input.vmag;
th0 = input.theta;
tflight = input.tflight;
nGrid = param.nGrid;

x0 = 0;
y0 = 0;
dx0 = v0*cos(th0);
dy0 = v0*sin(th0);
if dy0 < 0
    error('dy0 is less than 0!');
end

userFun = @(t, x) dyn_fun(t, x, param);
tSpan = [0, tflight];
x0 = [x0; y0; dx0; dy0];
solution = ode45(userFun, tSpan, x0);

sol.t = linspace(solution.x(1), solution.x(end), nGrid);
xx = deval(solution, sol.t);
sol.x = xx(1,:);
sol.y = xx(2,:);
sol.dx = xx(3,:);
sol.dy = xx(4,:);

end

function [cneq, ceq] = nl_cons(var_input, target, param)

x0  = 0;
y0  = 0;
dx0 = var_input(1);
dy0 = var_input(2);
T   = var_input(3);

nGrid = param.nGrid;
tSpan = linspace(0, T, nGrid);

userFun = @(t, x) dyn_fun(t, x, param);
x_0 = [x0; y0; dx0; dy0];
solution = ode113(userFun, tSpan, x_0);

t = linspace(solution.x(1), solution.x(end), nGrid);
xsol = deval(solution, t);

xfinal = xsol(1,end);
yfinal = xsol(2,end);

cneq = [];
ceq = [xfinal - target.x; yfinal - target.y];

end

function dx_vect = dyn_fun(~, x_vect, param)

dx_vect = zeros(size(x_vect));
dx_vect(1:2,:) = x_vect(3:4,:);

% Speed 
vx = x_vect(3,:);
vy = x_vect(4,:);
v = sqrt(vx.*vx + vy.*vy); 

% rho = 1.225;
% S= 1;
% Cd = 0.7;
% F_drag = @(v) 0.5*rho*norm(v)^2*S*Cd;

% Forces
% fx = -F_drag(vx);
fx = 0;
% fy = -param.g -F_drag(vy);
fy = -param.g;

% Derivative of velocity states
dx_vect(3,:) = fx; % unit point mass
dx_vect(4,:) = fy; % unit point mass

end

