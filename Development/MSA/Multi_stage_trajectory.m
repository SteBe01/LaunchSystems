clearvars; clc; close all;
clear ballisticMultiStagev2
clear ballisticMultiStage

% General parameters
stages.stage_1.n = 5.53;
stages.stage_2.n = 2.93;
stages.stage_1.Isp = 280;
stages.stage_2.Isp = 298;
stages.stage_1.Thrust = 252e3;
stages.stage_2.Thrust = 27.5e3;
stages.stage_1.m0 = 15e3;

stages.stage_1.d = 1.8;
stages.stage_2.d = 1.5;
stages.stage_1.Cd = 0.5;
stages.stage_2.Cd = 0.5;
stages.stage_1.Cl = 0.5;
stages.stage_2.Cl = 0.0;

params.g0 = 9.81;
params.Re = 6378000;
params.h0 = 11.9e3;
params.v0 = 400;
params.gamma0 = deg2rad(0);
params.turn_duration = 60;   % [s]
params.h_stage = 95e3;

% Pitch maneuver
params.t_turn = 5;                       % [m]    - Initial maneuver altitude
params.gamma_turn = deg2rad(30);        % [rad]   - Initial flight path angle

% MECO to stage separation wait time
stages.stage_1.t_wait = 5;

% Stage wait time before ignition
stages.stage_2.t_ign = 3;

% Compute necessary data
fn = fieldnames(stages);
for ii = 1:numel(fn)
    stages.(fn{ii}).m_prop = stages.(fn{ii}).m0 * (1 - 1/stages.(fn{ii}).n);
    stages.(fn{ii}).m_dot = stages.(fn{ii}).Thrust / (stages.(fn{ii}).Isp * params.g0);
    stages.(fn{ii}).t_burn = stages.(fn{ii}).m_prop / stages.(fn{ii}).m_dot;
    if ii < numel(fn)
        stages.(fn{ii+1}).m0 = stages.(fn{ii}).m0 - stages.(fn{ii}).m_prop;
    end
end

t_tot = stages.stage_1.t_burn + stages.stage_2.t_burn;
y0_s1 = [params.v0 params.gamma0 0 params.h0];
% y0_s1 = [0 params.h0 params.v0 0 params.gamma0];

options = odeset('RelTol',1e-8, 'MaxStep', 0.1, 'Events', @(t, y) stage_Separation(t, y, stages.stage_1));
[T1, Y1, te, ye, ie] = ode45(@(t, y) ballisticMultiStage(t, y, stages.stage_1, params, 1), [0:2*t_tot], y0_s1, options);

y0_s2 = Y1(end, :);
options = odeset('RelTol', 1e-8, 'MaxStep', 0.1, 'Events', @orbit_insertion);
[T2, Y2] = ode113(@(t, y) ballisticMultiStage(t, y, stages.stage_2, params, 2), [0:10*t_tot], y0_s2, options);

T = [T1; T2+T1(end)];
Y = [Y1; Y2];

% Retrieve data from ode
qdyn = zeros(length(T), 1);
gamma_dot = zeros(length(T), 1);
acc = zeros(length(T), 2);
for ii = 1:length(T)
    if ii <= length(T1)
        [~, parout] = ballisticMultiStage(T(ii), Y(ii, :), stages.stage_1, params, 1);
    else
        [~, parout] = ballisticMultiStage(T(ii), Y(ii, :), stages.stage_2, params, 2);
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

%% Presentation plots

% gamma_fig = figure("WindowStyle","normal"); hold on; grid on; title("$\mathbf{\gamma\ evolution\ in\ time}$", 'Interpreter','latex'); xlabel("$\mathbf{Time\ [s]}$", 'Interpreter','latex'); ylabel("$\mathbf{\gamma\ [deg]}$", 'Interpreter','latex');
% plot(T(1:end), rad2deg(Y(:, 2))); xlim([-5 70]); xline(param.t_turn, '--', "Manoeuvre start", "LabelVerticalAlignment","top", 'Interpreter','latex'); xline(param.t_turn+param.turn_duration, '--', "Manoeuvre ends", "LabelVerticalAlignment","top", 'Interpreter','latex');

% acc_fig = figure("WindowStyle","normal"); hold on; grid on; title("$\mathbf{\gamma\ evolution\ in\ time}$", 'Interpreter','latex'); xlabel("$\mathbf{Time\ [s]}$", 'Interpreter','latex'); ylabel("$\mathbf{\gamma\ [deg]}$", 'Interpreter','latex');
% plot(T(1:end), acc);

%% Plots for 3 dof with both body velocities (ballisticv2)

% set(0, 'DefaultLineLineWidth', 1.5)
% 
% boundary = 0;
% subplot(2,2,1), hold on, grid on, title("Downrange over time"), xlabel("Time [s]"), ylabel("Downrange [km]")
% plot(T, Y(1:end-boundary,1)/1e3)
% xline(T1(end), '--k', 'Staging');
% subplot(2,2,2), hold on, grid on, title("Altitude over time"), xlabel("Time [s]"), ylabel("Altitude [km]")
% plot(T, Y(1:end-boundary,2)/1e3)
% xline(T1(end), '--k', 'Staging');
% subplot(2,2,3), hold on, grid on, title("Axial velocity over time"), xlabel("Time [s]"), ylabel("Velocity [km/s]")
% plot(T, Y(1:end-boundary,3)/1e3)
% xline(T1(end), '--k', 'Staging');
% subplot(2,2,4), hold on, grid on, title("Trasversal velocity over time"), xlabel("Time [s]"), ylabel("Velocity [km/s]")
% plot(T, Y(1:end-boundary,4)/1e3)
% xline(T1(end), '--k', 'Staging');
% 
% figure; hold on; grid on; title("Dynamic pressure wrt altitude"), xlabel("Altitude [km]"), ylabel("Qdyn [kPa]")
% plot(Y(1:end, 2)/1e3, qdyn/1e3);
% xline(T1(end), '--k', 'Staging');
% 
% figure; hold on; grid on; title("Altitude wrt Downrange"), xlabel("Downrange [km]"), ylabel("Altitude [km]")
% plot(Y(1:end, 1)/1e3, Y(1:end, 2)/1e3);
% xline(T1(end), '--k', 'Staging');
% 
% figure; 
% subplot(2,1,1), hold on, grid on, title("Gamma over time"), xlabel("Time [s]"), ylabel("Gamma [rad]")
% plot(T, Y(1:end-boundary,5))
% xline(T1(end), '--k', 'Staging');
% subplot(2,1,2); hold on; grid on; title("Gamma dot over time"), xlabel("Time [s]"), ylabel("gamma dot [rad/s]")
% plot(T, gamma_dot);
% xline(T1(end), '--k', 'Staging');
% 
% figure;
% subplot(2,1,1); hold on; grid on; title("Axial acceleration over time"), xlabel("Time [s]"), ylabel("Acceleration [g]")
% plot(T, acc(:,1)/params.g0);
% xline(T1(end), '--k', 'Staging');
% subplot(2,1,2); hold on; grid on; title("Trasversal acceleration over time"), xlabel("Time [s]"), ylabel("Acceleration [g]")
% plot(T, acc(:,2)/params.g0);
% xline(T1(end), '--k', 'Staging');

%% Functions

function [dY, parout] = ballisticMultiStagev2(t,y, stage, params, current_stage)

    persistent turn_complete t_turn

    if isempty(turn_complete)
        turn_complete = false;
        t_turn = nan;
    end

    % Retrieve data from ode
    % downrange = y(1);
    altitude = y(2);
    vx_body = y(3);
    vy_body = y(4);
    gamma = y(5);

    if current_stage == 1
        t_wait = params.t_turn;
    else 
        t_wait = stage.t_ign;
    end

    % if current_stage == 1 && (t >= params.t_turn && t <= params.t_turn + params.turn_duration)
    if current_stage == 1 && t >= params.t_turn && (gamma - params.gamma_turn) < 1e-3 && ~turn_complete
    % if current_stage == 1 && t >= params.t_turn && gamma < params.gamma_turn && ~turn_complete
        delta = deg2rad(12);
    elseif current_stage == 1 && t >= params.t_turn && (gamma - params.gamma_turn) >= 1e-3 && ~turn_complete
    % elseif current_stage == 1 && t >= params.t_turn && gamma >= params.gamma_turn && ~turn_complete
        t_turn = t;
        turn_complete = true;
        delta = 0;
    else
        delta = 0;
    end
    velsNorm = norm([vx_body vy_body]);

    % Retrieve data used multiple times 
    t_burn = stage.t_burn;
    Re = params.Re;

    % [m^2] - Rocket surface area
    S = pi*(stage.d^2/4); 

    % [kg/m^3] - Density at current altitude
    rho = getAtmosphericData(altitude);

    if t <= t_wait
        m = stage.m0;
    elseif t <= t_burn + t_wait && t > t_wait
        m = stage.m0 - stage.m_dot * (t-t_wait);
    else
        m = stage.m0 - stage.m_dot * t_burn;
    end

    % [m/s^2] - Gravity acceleration taking into account altitude
    g = params.g0/((1+altitude/Re)^2);

    % [Pa] - Dynamic pressure
    qdyn = 0.5*rho*velsNorm^2;

    % [N] - Drag force acting on the rocket
    D = qdyn*S*stage.Cd;

    % [N] - Lift force acting on the rocket
    L = qdyn*S*stage.Cl;

    % [N] - Thrust force acting on the rocket
    if t <= t_burn + t_wait && t > t_wait
        T = stage.Thrust;
    else
        T = 0;
    end

    % Initialize derivative vector
    dY = zeros(5, 1);

    if current_stage == 1 && t < t_wait
        dY(1) = Re/(Re+altitude) * velsNorm;
        dY(2) = vy_body;
        dY(3) = -D/m;
        dY(4) = -g;
    else
        gammaDot = velsNorm*cos(gamma)/(Re+altitude) - g*cos(gamma)/velsNorm + T*sin(delta)/(m*velsNorm) + L/(m*velsNorm);
        dY(1) = Re/(Re+altitude) * (vx_body * cos(gamma) - vy_body * sin(gamma));
        dY(2) = vx_body * cos(gamma) + vy_body * sin(gamma);
        dY(3) = T/m*cos(delta) - D/m - g*sin(gamma); % - gammaDot*vy_body;
        dY(4) = T/m*sin(delta) - L/m - g*cos(gamma); % + gammaDot*vx_body;
        dY(5) = gammaDot;
    end

    % if current_stage == 1 && (t >= t_wait && t <= t_wait+params.turn_duration)
    %     dY(5) = (params.gamma_turn-0)*(pi/(2*params.turn_duration)*sin(pi*(t-params.t_turn)/params.turn_duration));
    % end

    % Prepare output struct for ode recall
    if nargout > 1
        parout.qdyn = qdyn;
        parout.gamma_dot = dY(5);
        parout.acc = reshape(dY(3:4), [1 2]);
        parout.t_turn = t_turn;
    end

end


function [dY, parout] = ballisticMultiStage(t,y, stage, params, current_stage)

    persistent turn_complete gamma_drop

    if isempty(turn_complete)
        turn_complete = false;
    end

    % Retrieve data from ode
    v = y(1);
    gamma = y(2);
    % x = y(3);
    h = y(4);

    if t < params.t_turn
        gamma_drop = gamma;
    end

    if current_stage == 1
        t_wait = params.t_turn;
    else 
        t_wait = stage.t_ign;
    end

    % Thrust vectoring

    % % if current_stage == 1 && (t >= params.t_turn && t <= params.t_turn + params.turn_duration)
    % if current_stage == 1 && t >= params.t_turn && (gamma - params.gamma_turn) < 1e-3 && ~turn_complete
    % % if current_stage == 1 && t >= params.t_turn && gamma < params.gamma_turn && ~turn_complete
    %     delta = -12;
    % elseif current_stage == 1 && t >= params.t_turn && (gamma - params.gamma_turn) >= 1e-3 && ~turn_complete
    % % elseif current_stage == 1 && t >= params.t_turn && gamma >= params.gamma_turn && ~turn_complete
    %     turn_complete = true;
    %     delta = 0;
    % else
        delta = 0;
    % end

    % Retrieve data used multiple times 
    t_burn = stage.t_burn;
    Re = params.Re;

    % [m^2] - Rocket surface area
    S = pi*(stage.d^2/4); 

    % [kg/m^3] - Density at current altitude
    rho = getAtmosphericData(h);

    if t <= t_wait
        m = stage.m0;
    elseif t <= t_burn + t_wait && t > t_wait
        m = stage.m0 - stage.m_dot * (t-t_wait);
    else
        m = stage.m0 - stage.m_dot * t_burn;
    end

    % [m/s^2] - Gravity acceleration taking into account altitude
    g = params.g0/((1+h/Re)^2);

    % [Pa] - Dynamic pressure
    qdyn = 0.5*rho*v^2;

    % [N] - Drag force acting on the rocket
    D = qdyn*S*stage.Cd;

    % [N] - Lift force acting on the rocket
    L = qdyn*S*stage.Cl;

    % [N] - Thrust force acting on the rocket
    if t <= t_burn + t_wait && t > t_wait
        T = stage.Thrust;
    else
        T = 0;
    end

    % Initialize derivative vector
    dY = zeros(4, 1);

    if current_stage == 1 && t < t_wait
        dY(1) = -D/m;
        dY(3) = Re/(Re+h) * v;
        dY(4) = -g*t;
    else
        dY(1) = T/m*cos(delta) - D/m - g*sin(gamma);
        dY(2) = v*cos(gamma)/(Re+h) - g*cos(gamma)/v + T*sin(delta)/(m*v) + L/(m*v);
        dY(3) = Re/(Re+h) * v * cos(gamma);
        dY(4) = v * sin(gamma);
    end

    if current_stage == 1 && (t >= params.t_turn && t <= params.t_turn+params.turn_duration)
        dY(2) = (params.gamma_turn-gamma_drop)*(pi/(2*params.turn_duration)*sin(pi*(t-params.t_turn)/params.turn_duration));
    end

    % Prepare output struct for ode recall
    if nargout > 1
        parout.qdyn = qdyn;
        parout.gamma_dot = dY(2);
        parout.acc = dY(1);
    end

end

function [value, isterminal, direction] = stage_Separation(t, ~, stage)
    value = t - (stage.t_burn + stage.t_wait);
    isterminal = 1;
    direction = 1;
end

function [value, isterminal, direction] = orbit_insertion(~, y)
    % value = y(4);% - 400e3;
    value = y(2) + deg2rad(2);
    isterminal = 1;
    direction = 0;
end
function [value, isterminal, direction] = orbit_insertionv2(~, y)
    value = y(2) - 400e3;
    isterminal = 1;
    direction = 0;
end