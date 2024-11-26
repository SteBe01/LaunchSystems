% TODO
% Implement the fuel derivative as done during classes

%% Simulator

clear, clc
close all
clear dyn

[stages, params, init] = loadMission();

% Compute masses and t_burn total
fn = fieldnames(stages);
for ii = 1:length(fn)
    stages.(fn{ii}).m_prop = stages.(fn{ii}).m0 * (1 - 1/stages.(fn{ii}).MR);
    stages.(fn{ii}).m_dot = stages.(fn{ii}).Thrust / (stages.(fn{ii}).Isp * params.g0);
    stages.(fn{ii}).t_burn_tot = stages.(fn{ii}).m_prop / stages.(fn{ii}).m_dot;
    if ii < numel(fn)
        stages.(fn{ii+1}).m0 = stages.(fn{ii}).m0 - stages.(fn{ii}).m_prop;
    end
end

% Trajectory propagation
y0_stg1 = [init.x0 init.z0 init.vx0 init.vz0 init.theta0 init.thetaDot0];

options_stg1 = odeset('RelTol',1e-8, 'MaxStep', 0.1, 'Events', @(t, y) stage_Separation(t, y, stages.stg1));
options_stg2 = odeset('RelTol', 1e-8, 'MaxStep', 0.1, 'Events', @(t, y) orbit_insertion(t, y));

tspan = 0:1/params.u_freq:1e4;

T1 = zeros(length(tspan), 1);
Y1 = zeros(length(tspan), 6);
delta_vec = zeros(length(tspan), 1);

Y1(1, :) = y0_stg1;
idx = 2;

alpha = 0;
theta = 0;
thetaDot = 0;


% First stage simulation
for ii = 2:length(tspan)-1

    if Y1(idx-1, 2) < 50e3
        angle = 90;
    else
        angle = 45;
    end
    err = theta - deg2rad(angle);

    delta = -params.k1*err - params.k2*thetaDot - params.k3*alpha;
    if abs(delta) > stages.stg1.deltaMax 
        delta = stages.stg1.deltaMax*sign(delta);
    end
    
    [tvect, yvect, tevent, yevent,~] = ode113(@(t,y) dyn(t, y, stages.stg1, params, 1, delta), [tspan(ii) tspan(ii+1)], Y1(idx-1,:), options_stg1);
    parout = recallOdeFcn(tvect, yvect, stages.stg1, params, 1, delta);
    T1(idx:idx+length(tvect)-1) = tvect;
    Y1(idx:idx+length(tvect)-1, :) = yvect;
    delta_vec(idx:idx+length(tvect)-1) = repmat(delta, [length(tvect) 1]);
    idx = idx+length(tvect);

    alpha = parout.alpha(end);
    theta = yvect(end, 5);
    thetaDot = yvect(end, 6);

    % Exit the loop if we have stage separation
    if ~isempty(tevent)
        break
    end
end
T1(idx:end, :) = [];
Y1(idx:end, :) = [];
delta_vec(idx:end, :) = [];
clear dyn

% [T1, Y1] = ode113(@(t, y) dyn(t, y, stages.stg1, params, 1), 0:1e4, y0_stg1, options_stg1);
y0_s2 = Y1(end, :);
% clear dyn
[T2, Y2] = ode113(@(t, y) dyn(t, y, stages.stg2, params, 2, 0), 0:1e4, y0_s2, options_stg2);
clear dyn

T = [T1; T2+T1(end)];
Y = [Y1; Y2];
delta_vec = [delta_vec; zeros(length(T2), 1)];

% Retrieve data from ode
qdyn = zeros(length(T), 1);
acc = zeros(length(T), 2);
alpha = zeros(length(T), 1);
moment = zeros(length(T), 1);
for ii = 1:length(T)
    if ii <= length(T1)
        [~, parout] = dyn(T(ii), Y(ii, :), stages.stg1, params, 1, delta_vec(ii));
    else
        [~, parout] = dyn(T(ii), Y(ii, :), stages.stg2, params, 2, 0);
    end
    qdyn(ii) = parout.qdyn;
    acc(ii,:) = parout.acc;
    % if isfield(parout, "t_turn") && ~isnan(parout.t_turn)
    %     t_turn = parout.t_turn;
    % end
    alpha(ii) = parout.alpha;
    moment(ii) = parout.moment;
end

downrange = params.Re./(params.Re+Y(:,1)) .* Y(:,1);


%% Plots for simple 3 dof

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
plot(T, unwrap(rad2deg(alpha)))
xline(T1(end), '--k', 'Staging')
subplot(4,1,4), hold on, grid on, title("$\delta$ evolution", 'Interpreter', 'latex'), xlabel("Time [s]"), ylabel("Delta [deg]")
plot(T, unwrap(rad2deg(delta_vec)))
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


%% Functions

function [value, isterminal, direction] = stage_Separation(t, ~, stage)
    value = t - (stage.t_burn_tot + stage.t_wait);
    isterminal = 1;
    direction = 1;
end

function [value, isterminal, direction] = orbit_insertion(~, y)
    value = y(4);
    isterminal = 1;
    direction = 0;
end

