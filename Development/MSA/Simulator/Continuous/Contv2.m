clear, clc
close all
clear dyn stage_Separation orbit_insertion

addpath(genpath("Functions"))
addpath(genpath("Functions_events"))

[stages, params, init] = loadMission();

stop = false;

ii = 1;
T1 = [];
Y1 = [];
Y0 = [init.x0 init.z0 init.vx0 init.vz0 init.theta0 init.thetaDot0 init.m_prop];
t0 = 0;
current_altitude = params.pitch.initial_altitude;
params.dispStat = false;

angle = deg2rad(46);
Dangle = 0;

while ~stop 
    angle = angle+Dangle;
    params.pitch.first_angle = angle;
    params.pitch.initial_altitude = current_altitude;

    [T, Y, ie] = simulate(stages, params, t0, Y0, 1);
    
    r = norm(Y(end,1:2))-params.Re;
    if r < 399.5e3 && ie ~= 2
        Dangle = deg2rad(abs(r-399.5e3)/(25e3/4.5));
    elseif r > 400.5e3
        Dangle = -deg2rad(abs(r-400.5e3)/(25e3/4.5));
    elseif ie==2
        Dangle = deg2rad(5* rand());
    else
        Dangle = 0;
        stop = 1;
    end
    if abs(Dangle) > deg2rad(5)
        Dangle = deg2rad(5) * sign(Dangle);
    end

    disp("Current Error: " + num2str(r/1e3 - 400) + ", with angle: " + num2str(rad2deg(angle)) + " and velocities: [" + num2str(Y(end,3)/1e3) + ", " + num2str(Y(end,4)/1e3) + "]");
end

fprintf("Reached the following condition: angle = %.3f, error: %.3f, velocities: [%.3f, %.3f]\n", rad2deg(angle), r/1e3 - 400, Y(end,3)/1e3, Y(end,4)/1e3);

%%

params.dispStat = true;
params.pitch.first_angle = angle;
params.pitch.initial_altitude = current_altitude;

[T, Y, idxStage, parout] = run_simulator(stages, params, init, 1);
plotData(T, Y, params, parout, idxStage);


%% Auxiliary function

function [T, Y, ie] = simulate(stages, params, t0, Y0, full_flight)

    clear dyn stage_Separation orbit_revolution

    t_max = 1e3;
    
    %% First stage simulation

    options_stg1 = odeset('RelTol',1e-8, 'MaxStep', 0.1, 'Events', @(t, y) stage_Separation(t, y, stages.stg1, params));

    [T1, Y1, ~, ~, ie] = ode45(@(t,y) dyn(t, y, stages.stg1, params, 1), [t0 t_max], Y0, options_stg1);
    clear dyn

    T = T1;
    Y = Y1;

    %% Second stage simulation
    if full_flight && ie~=2
        options_stg2 = odeset('RelTol', 1e-8, 'MaxStep', 0.1, 'Events', @(t, y) orbit_insertion(t, y));
        [T2, Y2] = ode45(@(t,y) dyn(t, y, stages.stg2, params, 2), [0 t_max], [Y1(end,1:end-1) stages.stg2.m_prop], options_stg2);
        clear dyn

        T = [T; T2+T1(end)];
        Y = [Y; Y2];
    end
end