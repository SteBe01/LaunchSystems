%% Simulator

clear, clc
close all
clear dyn

[stages, params, init] = loadMission();

% Trajectory propagation
y0_stg1 = [init.x0 init.z0 init.vx0 init.vz0 init.theta0 init.thetaDot0];

options_stg1 = odeset('RelTol',1e-8, 'MaxStep', 0.1, 'Events', @(t, y) stage_Separation(t, y, stages.stg1));
options_stg2 = odeset('RelTol', 1e-8, 'MaxStep', 0.1, 'Events', @(t, y) orbit_insertion(t, y));

t_max = 1e4;
% First stage simulation
[T1, Y1] = ode113(@(t,y) dyn(t, y, stages.stg1, params, 1), [0 t_max], y0_stg1, options_stg1);
clear dyn

% Second stage simulation
[T2, Y2] = ode113(@(t,y) dyn(t, y, stages.stg2, params, 2), [0 t_max], Y1(end,:), options_stg2);
clear dyn


%% Retrieve data from ode

T = [T1; T2+T1(end)];
Y = [Y1; Y2];

parout_stg1 = recallOdeFcn(T1, Y1, stages.stg1, params, 1);
parout_stg2 = recallOdeFcn(T2, Y2, stages.stg2, params, 2);

qdyn = [parout_stg1.qdyn; parout_stg2.qdyn];
acc = [parout_stg1.acc; parout_stg2.acc];
alpha = [parout_stg1.alpha; parout_stg2.alpha];
moment = [parout_stg1.moment; parout_stg2.moment];
dv_drag_vec = [parout_stg1.dv_drag_vec; parout_stg2.dv_drag_vec];
dv_grav_vec = [parout_stg1.dv_grav_vec; parout_stg2.dv_grav_vec];
delta_vec = [parout_stg1.delta; parout_stg2.delta];

dv_drag_s1 = cumtrapz(parout_stg1.dv_drag_vec, T1);
dv_grav_s1 = cumtrapz(parout_stg1.dv_grav_vec, T1);
dv_thrust_s1 = stages.stg1.Isp*9.81*log(stages.stg1.m0/parout_stg1.m(end));
dv_s1 = dv_drag_s1(end) + dv_grav_s1(end) + dv_thrust_s1;

dv_drag_s2 = cumtrapz(parout_stg2.dv_drag_vec, T2);
dv_grav_s2 = cumtrapz(parout_stg2.dv_grav_vec, T2);
dv_thrust_s2 = stages.stg2.Isp*9.81*log(stages.stg2.m0/parout_stg2.m(end));
dv_s2 = dv_drag_s2(end) + dv_grav_s2(end) + dv_thrust_s2;

g_vec = params.g0./((1+Y1(:,2)/params.Re).^2);
downrange = params.Re./(params.Re+Y(:,1)) .* Y(:,1);


%% Plots 

set(0, 'DefaultLineLineWidth', 1.5)

boundary = 0;
subplot(2,2,1), hold on, grid on, title("Downrange over time"), xlabel("Time [s]"), ylabel("Downrange [km]")
plot(T, downrange/1e3)
xline(T1(end), '--k', 'Staging')
subplot(2,2,2), hold on, grid on, title("Vertical position over time"), xlabel("Time [s]"), ylabel("Altitude [km]")
plot(T, Y(1:end-boundary,2)/1e3)
xline(T1(end), '--k', 'Staging')
subplot(2,2,3), hold on, grid on, title("Horizontal velocity over time"), xlabel("Time [s]"), ylabel("Velocity [km/s]")
plot(T, Y(1:end-boundary,3)/1e3)
xline(T1(end), '--k', 'Staging')
subplot(2,2,4), hold on, grid on, title("Vertical velocity over time"), xlabel("Time [s]"), ylabel("Velocity [km/s]")
plot(T, Y(1:end-boundary,4)/1e3)
xline(T1(end), '--k', 'Staging')

figure, hold on, grid on, title("Dynamic pressure wrt altitude"), xlabel("Altitude [km]"), ylabel("Qdyn [kPa]")
plot(Y(1:end, 2)/1e3, qdyn/1e3);
xline(Y1(end, 1)/1e3, '--k', 'Staging')

figure, hold on, grid on, title("Altitude wrt Downrange"), xlabel("Downrange [km]"), ylabel("Altitude [km]")
plot(downrange/1e3, Y(1:end, 2)/1e3)
xline(downrange(length(T1))/1e3, '--k', 'Staging')

figure
subplot(4,1,1), hold on, grid on, title("Theta over time"), xlabel("Time [s]"), ylabel("Theta [deg]")
plot(T, rad2deg(Y(:, 5)))
xline(T1(end), '--k', 'Staging')
subplot(4,1,2), hold on, grid on, title("Theta dot over time"), xlabel("Time [s]"), ylabel("Theta dot [deg/s]")
plot(T, rad2deg(Y(:, 6)))
xline(T1(end), '--k', 'Staging')
subplot(4,1,3), hold on, grid on, title("$\alpha$ evolution", 'Interpreter', 'latex'), xlabel("Time [s]"), ylabel("Alpha [deg]")
plot(T, rad2deg(alpha))
xline(T1(end), '--k', 'Staging')
subplot(4,1,4), hold on, grid on, title("$\delta$ evolution", 'Interpreter', 'latex'), xlabel("Time [s]"), ylabel("Delta [deg]")
plot(T, rad2deg(delta_vec))
xline(T1(end), '--k', 'Staging')

figure
subplot(2,2,1), hold on, grid on, title("Axial acceleration over time"), xlabel("Time [s]"), ylabel("Acceleration [g]")
plot(T, acc(:,2)/params.g0)
xline(T1(end), '--k', 'Staging')
subplot(2,2,2), hold on, grid on, title("Normal acceleration over time"), xlabel("Time [s]"), ylabel("Acceleration [g]")
plot(T, acc(:,1)/params.g0)
plot(T, -1./((1+Y(:,2)./params.Re).^2), '--')
xline(T1(end), '--k', 'Staging')
subplot(2,2,3), hold on, grid on, title("Acceleration norm over time"), xlabel("Time [s]"), ylabel("Acceleration norm [g]")
plot(T, vecnorm(acc, 2,2)./params.g0)
plot(T, 1./((1+Y(:,2)./params.Re).^2), '--')
xline(T1(end), '--k', 'Staging')
subplot(2,2,4), hold on, grid on, title("Moments over time"), xlabel("Time [s]"), ylabel("Moment [Nm]")
plot(T, moment)
xline(T1(end), '--k', 'Staging')


%% Event functions

function [value, isterminal, direction] = stage_Separation(t, y, stage)
    value = t - (stage.t_burn_tot + stage.t_wait + 1);
    % value = y(2);
    isterminal = 1;
    direction = 0;
end

function [value, isterminal, direction] = orbit_insertion(~, y)
    value = y(4);
    isterminal = 1;
    direction = 0;
end

