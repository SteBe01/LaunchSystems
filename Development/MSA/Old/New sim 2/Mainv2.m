% TODO
% Implement the fuel derivative as done during classes

%% Simulator

clear, clc
close all

[stages, params] = loadMission();

% Compute masses and t_burn total
fn = fieldnames(stages);
for ii = 1:length(fn)
    stages.(fn{ii}).m_prop = stages.(fn{ii}).m0 * (1 - 1/stages.(fn{ii}).n);
    stages.(fn{ii}).m_dot = stages.(fn{ii}).Thrust / (stages.(fn{ii}).Isp * params.g0);
    stages.(fn{ii}).t_burn_tot = stages.(fn{ii}).m_prop / stages.(fn{ii}).m_dot;
    if ii < numel(fn)
        stages.(fn{ii+1}).m0 = stages.(fn{ii}).m0 - stages.(fn{ii}).m_prop;
    end
end
t_tot = stages.stg1.t_burn_tot + stages.stg2.t_burn_tot;

% Trajectory propagation
y0_stg1 = [0 params.h0 params.v0 0 params.gamma0];

options_stg1 = odeset('RelTol',1e-8, 'MaxStep', 0.1, 'Events', @(t, y) stage_Separation(t, y, stages.stg1));
options_stg2 = odeset('RelTol', 1e-8, 'MaxStep', 0.1, 'Events', @(t, y) orbit_insertion(t, y));

[T1, Y1] = ode113(@(t, y) rocket_dynamics(t, y, stages.stg1, params, 1), [0 1e5], y0_stg1, options_stg1);
y0_s2 = Y1(end, :);
[T2, Y2] = ode113(@(t, y) rocket_dynamics(t, y, stages.stg2, params, 2), [0 1e3], y0_s2, options_stg2);
clear rocket_dynamics

T = [T1; T2+T1(end)];
Y = [Y1; Y2];

% Retrieve data from ode
qdyn = zeros(length(T), 1);
gammaDot = zeros(length(T), 1);
acc = zeros(length(T), 2);
Mach = zeros(length(T), 1);
AoA = zeros(length(T), 1);
dcm = zeros(2,2, length(T));
for ii = 1:length(T)
    if ii <= length(T1)
        [~, parout] = rocket_dynamics(T(ii), Y(ii, :), stages.stg1, params, 1);
    else
        [~, parout] = rocket_dynamics(T(ii), Y(ii, :), stages.stg2, params, 2);
    end
    qdyn(ii) = parout.qdyn;
    gammaDot(ii) = parout.gammaDot;
    acc(ii,:) = parout.acc;
    Mach(ii,:) = parout.Mach;
    AoA(ii,:) = parout.AoA;
    if isfield(parout, "t_turn") && ~isnan(parout.t_turn)
        t_turn = parout.t_turn;
    end
    % dcm(:,:,ii) = parout.dcm;
end

vel_body = cumtrapz(T, acc);


%% Plots for simple 3 dof

set(0, 'DefaultLineLineWidth', 1.5)

boundary = 0;
subplot(2,2,1), hold on, grid on, title("Downrange over time"), xlabel("Time [s]"), ylabel("Downrange [km]")
plot(T, Y(1:end-boundary,1)/1e3)
xline(T1(end), '--k', 'Staging')
subplot(2,2,2), hold on, grid on, title("Altitude over time"), xlabel("Time [s]"), ylabel("Altitude [km]")
plot(T, Y(1:end-boundary,2)/1e3)
xline(T1(end), '--k', 'Staging')
subplot(2,2,3), hold on, grid on, title("Horizontal velocity over time"), xlabel("Time [s]"), ylabel("Velocity [km/s]")
plot(T, Y(1:end-boundary,3)/1e3)
xline(T1(end), '--k', 'Staging')
subplot(2,2,4), hold on, grid on, title("Vertical velocity over time"), xlabel("Time [s]"), ylabel("Velocity [deg]")
plot(T, Y(1:end-boundary,4)/1e3)
xline(T1(end), '--k', 'Staging')

figure; hold on; grid on; title("Dynamic pressure wrt altitude"), xlabel("Altitude [km]"), ylabel("Qdyn [kPa]")
plot(Y(1:end, 2)/1e3, qdyn/1e3)
xline(Y1(end, 2)/1e3, '--k', 'Staging')

figure; hold on; grid on; title("Altitude wrt Downrange"), xlabel("Downrange [km]"), ylabel("Altitude [km]")
plot(Y(1:end, 1)/1e3, Y(1:end, 2)/1e3)
xline(Y1(end,1)/1e3, '--k', 'Staging')

figure; 
subplot(2,1,1); hold on; grid on; title("Gamma over time"), xlabel("Time [s]"), ylabel("gamma [rad]")
plot(T, Y(:,5))
xline(T1(end), '--k', 'Staging')
subplot(2,1,2); hold on; grid on; title("Gamma dot over time"), xlabel("Time [s]"), ylabel("gamma\_dot [rad/s]")
plot(T, gammaDot)
xline(T1(end), '--k', 'Staging')

figure
subplot(2,1,1); hold on; grid on; title("Axial acceleration over time"), xlabel("Time [s]"), ylabel("Acceleration [g]")
plot(T, acc(:,1)/params.g0)
xline(T1(end), '--k', 'Staging')
subplot(2,1,2); hold on; grid on; title("Normal acceleration over time"), xlabel("Time [s]"), ylabel("Acceleration [g]")
plot(T, acc(:,2)/params.g0)
xline(T1(end), '--k', 'Staging')

figure
subplot(2,2,1), hold on, grid on, title("Horizontal velocity over time"), xlabel("Time [s]"), ylabel("Velocity [km/s]")
plot(T, Y(1:end-boundary,3)/1e3)
xline(T1(end), '--k', 'Staging')
subplot(2,2,2), hold on, grid on, title("Vertical velocity over time"), xlabel("Time [s]"), ylabel("Velocity [deg]")
plot(T, Y(1:end-boundary,4)/1e3)
xline(T1(end), '--k', 'Staging')
subplot(2,2,3), hold on, grid on, title("Axial velocity over time"), xlabel("Time [s]"), ylabel("Velocity [km/s]")
plot(T, vel_body(:,1)/1e3)
xline(T1(end), '--k', 'Staging')
subplot(2,2,4), hold on, grid on, title("Normal velocity over time"), xlabel("Time [s]"), ylabel("Velocity [deg]")
plot(T, vel_body(:,2)/1e3)
xline(T1(end), '--k', 'Staging')


%% Post-processing

% vel_ned = zeros(length(T), 2);
% for ii = 1:length(T)
%     vel_ned(ii,:) = dcm(:,:,ii) * vel_body(ii,:)' + [params.v0; 0];
% end
% 
% figure;
% subplot(2,1,1); hold on; grid on;
% plot(T, vel_ned(:,1)/1e3, 'DisplayName', 'recomputed');
% plot(T, Y(:, 3)/1e3, 'DisplayName', 'ODE');
% xline(T1(end), '--k', 'Staging', 'HandleVisibility','off');
% legend();
% subplot(2,1,2); hold on; grid on;
% plot(T, vel_ned(:,2)/1e3, 'DisplayName', 'recomputed');
% plot(T, Y(:, 4)/1e3, 'DisplayName', 'ODE');
% xline(T1(end), '--k', 'Staging', 'HandleVisibility', 'off');
% legend()


%% Event Functions

function [value, isterminal, direction] = stage_Separation(t, ~, stage)
    value = t - (stage.t_burn_tot + stage.t_wait);
    isterminal = 1;
    direction = 1;
end

function [value, isterminal, direction] = orbit_insertion(~, y)
    % value = y(2);% - 400e3;
    % value = y(5) + deg2rad(2);
    value = y(4);
    isterminal = 1;
    direction = 0;
end

