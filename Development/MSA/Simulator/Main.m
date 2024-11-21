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
y0_stg1 = [params.v0 params.gamma0 0 params.h0];

options_stg1 = odeset('RelTol',1e-8, 'MaxStep', 0.1, 'Events', @(t, y) stage_Separation(t, y, stages.stg1));
options_stg2 = odeset('RelTol', 1e-8, 'MaxStep', 0.1, 'Events', @(t, y) orbit_insertion(t, y));

[T1, Y1] = ode113(@(t, y) dyn(t, y, stages.stg1, params, 1), 0:2*t_tot, y0_stg1, options_stg1);
y0_s2 = Y1(end, :);
[T2, Y2] = ode113(@(t, y) dyn(t, y, stages.stg2, params, 2), 0:10*t_tot, y0_s2, options_stg2);

T = [T1; T2+T1(end)];
Y = [Y1; Y2];

% Retrieve data from ode
qdyn = zeros(length(T), 1);
gamma_dot = zeros(length(T), 1);
acc = zeros(length(T), 2);
for ii = 1:length(T)
    if ii <= length(T1)
        [~, parout] = dyn(T(ii), Y(ii, :), stages.stg1, params, 1);
    else
        [~, parout] = dyn(T(ii), Y(ii, :), stages.stg2, params, 2);
    end
    qdyn(ii) = parout.qdyn;
    gamma_dot(ii) = parout.gamma_dot;
    acc(ii,:) = parout.acc;
    if isfield(parout, "t_turn") && ~isnan(parout.t_turn)
        t_turn = parout.t_turn;
    end
end


%% Plots for simple 3 dof

set(0, 'DefaultLineLineWidth', 1.5)

boundary = 0;
subplot(2,2,1), hold on, grid on, title("Velocity over time"), xlabel("Time [s]"), ylabel("Velocity [km/s]")
plot(T, Y(1:end-boundary,1)/1e3)
xline(T1(end), '--k', 'Staging');
subplot(2,2,2), hold on, grid on, title("Gamma over time"), xlabel("Time [s]"), ylabel("Gamma [deg]")
plot(T, rad2deg(Y(1:end-boundary,2)))
xline(T1(end), '--k', 'Staging');
subplot(2,2,3), hold on, grid on, title("Downrange over time"), xlabel("Time [s]"), ylabel("Downrange [km]")
plot(T, Y(1:end-boundary,3)/1e3)
xline(T1(end), '--k', 'Staging');
subplot(2,2,4), hold on, grid on, title("Altitude over time"), xlabel("Time [s]"), ylabel("Altitude [km]")
plot(T, Y(1:end-boundary,4)/1e3)
xline(T1(end), '--k', 'Staging');

figure; hold on; grid on; title("Dynamic pressure wrt altitude"), xlabel("Altitude [km]"), ylabel("Qdyn [kPa]")
plot(Y(1:end, 4)/1e3, qdyn/1e3);
xline(T1(end), '--k', 'Staging');

figure; hold on; grid on; title("Altitude wrt Downrange"), xlabel("Downrange [km]"), ylabel("Altitude [km]")
plot(Y(1:end, 3)/1e3, Y(1:end, 4)/1e3)
xline(T1(end), '--k', 'Staging');

figure; hold on; grid on; title("Gamma dot over time"), xlabel("Time [s]"), ylabel("gamma_dot [rad/s]")
plot(T, gamma_dot);
xline(T1(end), '--k', 'Staging');

figure; hold on; grid on; title("Longitudinal acceleration over time"), xlabel("Time [s]"), ylabel("Acceleration [g]")
plot(T, acc(:,1)/params.g0);
xline(T1(end), '--k', 'Staging');


%% Functions

function [value, isterminal, direction] = stage_Separation(t, ~, stage)
    value = t - (stage.t_burn_tot + stage.t_wait);
    isterminal = 1;
    direction = 1;
end

function [value, isterminal, direction] = orbit_insertion(~, y)
    % value = y(4);% - 400e3;
    value = y(2) + deg2rad(2);
    isterminal = 1;
    direction = 0;
end


