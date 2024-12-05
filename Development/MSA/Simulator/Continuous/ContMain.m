%% Simulator

clear, clc
close all
clear dyn stage_Separation orbit_insertion

addpath(genpath("Functions"))
addpath(genpath("Functions_events"))

[stages, params, init] = loadMission();

angle = deg2rad(45);
t_vect = linspace(0,300,1000);

stop = 0;
ii = 1;
T1 = [];
Y1 = [];
y0_stg1 = [init.x0 init.z0 init.vx0 init.vz0 init.theta0 init.thetaDot0 init.m_prop];

while ~stop && ii < length(t_vect)-1
    options_stg1 = odeset('RelTol',1e-8, 'MaxStep', 0.1, 'Events', @(t, y) stage_Separation(t, y, stages.stg1, params, angle));
    [T1_temp, Y1_temp, ~, ~, ie] = ode45(@(t,y) dyn(t, y, stages.stg1, params, 1, angle), [t_vect(ii) t_vect(ii+1)], y0_stg1, options_stg1);
    % clear dyn
    T1 = [T1; T1_temp];
    Y1 = [Y1; Y1_temp];
    y0_stg1 = Y1(end,:);

    if (norm(Y1(end,1:2)) - 400e3) < 1000
        stop = 1;
    end
    if ie == 2
        stop = 1;
    end
    
    ii = ii + 1;
end


    % options_stg2 = odeset('RelTol', 1e-8, 'MaxStep', 0.1, 'Events', @(t, y) orbit_insertion(t, y)); %@(t,y) orbit_revolution(t, y, params));
    % [T2, Y2] = ode45(@(t,y) dyn(t, y, stages.stg2, params, 2, angle), [0 t_max], [Y1(end,1:end-1) stages.stg2.m_prop], options_stg2);
    % % clear dyn
    % 
    % T = [T; T2+T1(end)];
    % Y = [Y; Y2];

