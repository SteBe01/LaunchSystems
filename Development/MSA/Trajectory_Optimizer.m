clear, clc
close all

% General parameters
n = 7;                      % [-]         - Mass ratio m0/mf
I_sp = 309;                 % [s]         - Specific impulse
gamma0 = deg2rad(0);    % [rad]       - Initial pitch angle
param.m0 = 68e3;           % [kg]        - Initial mass 
param.g0 = 9.81;            % [m/s^2]     - Earth's gravitational parameter
param.Re = 6378000;         % [m]         - Earth's radius
param.Thrust = 933.91e3;    % [N]         - Max Thrust
param.d = 1.2;                % [m]         - Rocket's diameter
param.TW = 1.4;             % [-]         - Thrust over Weight ratio
param.Cd = 0.5;             % [-]         - Drag coefficient
param.rho0 = 1.225;         % [kg/m^3]    - Air density at sea level
param.H0 = 7500;            % [m]         - Air density reference altitude
param.gamma0 = gamma0;
param.turn_duration = 60;   % [s]


% Pitch maneuver
param.t_turn = 5;                       % [m]     - Initial maneuver altitude
param.gamma_turn = deg2rad(89.85);        % [rad]   - Initial flight path angle

% Compute necessary data
m_prop = param.m0*(1-1/n);        % [Kg]    - Propellant mass
param.m_dot = param.Thrust /(I_sp*param.g0);  % [Kg/s]  - Mass flow rate
param.t_burn = m_prop / param.m_dot;    % [s]     - Burning time

% ODE Initial state - V0 ~= 0 to avoid integration problems
y0 = [30 gamma0 0 12e3];
% options = odeset('RelTol',1e-8,'Events',@(x,y) event(x,y,1));
options = odeset('RelTol',1e-8);
[T,Y] = ode113(@(t, y) ballistic(t, y, param),[0:0.01:2*param.t_burn],y0,options);

% Retrieve data from ode
qdyn = zeros(length(T), 1);
gamma_dot = zeros(length(T), 1);
for ii = 1:length(T)
    [~, parout] = ballistic(T(ii), Y(ii, :), param);
    qdyn(ii) = parout.qdyn;
    gamma_dot(ii) = parout.gamma_dot;
    acc(ii) = parout.acc;
end

%% Plots

set(0, 'DefaultLineLineWidth', 1.5)

boundary = 0;
subplot(2,2,1), hold on, grid on, title("Velocity over time"), xlabel("Time [s]"), ylabel("Velocity [km/s]")
plot(T(1:end-boundary), Y(1:end-boundary,1)/1e3)
subplot(2,2,2), hold on, grid on, title("Gamma over time"), xlabel("Time [s]"), ylabel("Gamma [rad]")
plot(T(1:end-boundary), Y(1:end-boundary,2))
subplot(2,2,3), hold on, grid on, title("Downrange over time"), xlabel("Time [s]"), ylabel("Downrange [km]")
plot(T(1:end-boundary), Y(1:end-boundary,3)/1e3)
subplot(2,2,4), hold on, grid on, title("Altitude over time"), xlabel("Time [s]"), ylabel("Altitude [km]")
plot(T(1:end-boundary), Y(1:end-boundary,4)/1e3)

figure; hold on; grid on; title("Dynamic pressure wrt altitude"), xlabel("Altitude [km]"), ylabel("Qdyn [kPa]")
plot(Y(1:end, 4)/1e3, qdyn/1e3);

figure; hold on; grid on; title("Altitude wrt Downrange"), xlabel("Downrange [km]"), ylabel("Altitude [km]")
plot(Y(1:end, 3)/1e3, Y(1:end, 4)/1e3)

figure; hold on; grid on; title("Gamma dot over time"), xlabel("Time [s]"), ylabel("gamma_dot [rad/s]")
plot(T(1:end-boundary), gamma_dot);

%% Presentation plots

gamma_fig = figure("WindowStyle","normal"); hold on; grid on; title("$\mathbf{\gamma\ evolution\ in\ time}$", 'Interpreter','latex'); xlabel("$\mathbf{Time\ [s]}$", 'Interpreter','latex'); ylabel("$\mathbf{\gamma\ [deg]}$", 'Interpreter','latex');
plot(T(1:end), rad2deg(Y(:, 2))); xlim([-5 70]); xline(param.t_turn, '--', "Manoeuvre start", "LabelVerticalAlignment","top", 'Interpreter','latex'); xline(param.t_turn+param.turn_duration, '--', "Manoeuvre ends", "LabelVerticalAlignment","top", 'Interpreter','latex');

% acc_fig = figure("WindowStyle","normal"); hold on; grid on; title("$\mathbf{\gamma\ evolution\ in\ time}$", 'Interpreter','latex'); xlabel("$\mathbf{Time\ [s]}$", 'Interpreter','latex'); ylabel("$\mathbf{\gamma\ [deg]}$", 'Interpreter','latex');
% plot(T(1:end), acc);

%% Functions

% function [value, isterminal, direction] = event(~,xx,isTerminal)
%     value = xx(4);
%     isterminal = isTerminal;
%     direction = 0;
% end

function [dY, parout] = ballistic(t,y, param)

    persistent gamma_drop

    % Retrieve data from ode
    v = y(1);
    gamma = y(2);
    x = y(3);
    h = y(4);

    % Retrieve data used multiple times from param structure
    t_burn = param.t_burn;
    Re = param.Re;
    gamma0 = param.gamma0;

    % [m^2] - Rocket surface area
    S = pi*(param.d^2/4); 

    % [kg/m^3] - Density at current altitude
    rho = param.rho0*exp(-h/param.H0);

    % [kg] - Current rocket mass
    m = param.m0 - ((param.m_dot*t) * (t <= t_burn) + (param.m_dot*t_burn) * (t > t_burn));

    % [m/s^2] - Gravity acceleration taking into account altitude
    g = param.g0/((1+h/Re)^2);

    % [Pa] - Dynamic pressure
    qdyn = 0.5*rho*v^2;

    % [N] - Drag force acting on the rocket
    D = qdyn*S*param.Cd;

    % [N] - Thrust force acting on the rocket
    T = param.Thrust * (t < t_burn+param.t_turn) * (t > param.t_turn);

    % Initialize derivative vector
    dY = zeros(4, 1);

    % if t < param.t_turn
    %     dY(1) = T/m - D/m - g*sin(gamma);
    %     dY(2) = -1/v * (g - v^2/(Re+h)) * cos(gamma);
    %     dY(3) = Re/(Re+h) * v * cos(gamma);
    %     dY(4) = v * sin(gamma);
    % else
        dY(1) = T/m - D/m - g*sin(gamma);
        dY(2) = -1/v * (g - v^2/(Re+h)) * cos(gamma) * (t > param.t_turn);
        dY(3) = Re/(Re+h) * v * cos(gamma);
        dY(4) = v * sin(gamma);
    % end

    if t < param.t_turn
        gamma_drop = gamma;
    end
    if t >= param.t_turn && t <= param.t_turn+param.turn_duration
        dY(2) = (param.gamma_turn-gamma_drop)*(pi/(2*param.turn_duration)*sin(pi*(t-param.t_turn)/param.turn_duration));
    end

    % Prepare output struct for ode recall
    if nargout > 1
        parout.qdyn = qdyn;
        parout.gamma_dot = dY(2);
        parout.acc = dY(1);
    end

end

function [dY, parout] = ballistic_v2(t,y, param)

    % Retrieve data from ode
    vx = y(1);
    vy = y(2);
    gamma = y(3);
    x = y(4);
    h = y(5);

    % Retrieve data used multiple times from param structure
    t_burn = param.t_burn;
    Re = param.Re;

    % [m^2] - Rocket surface area
    S = pi*(param.d^2/4); 

    % [kg/m^3] - Density at current altitude
    [~, a, ~, rho] = atmosisa(h, "extended","on");

    % [kg] - Current rocket mass
    m = param.m0 - ((param.m_dot*t) * (t <= t_burn) + (param.m_dot*t_burn) * (t > t_burn));

    % [m/s^2] - Gravity acceleration taking into account altitude
    g = param.g0/((1+h/Re)^2);

    % [Pa] - Dynamic pressure
    qdyn = 0.5*rho*norm([vx vy])^2;

    % [N] - Drag force acting on the rocket
    D = qdyn*S*param.Cd;

    % [N] - Thrust force acting on the rocket
    T = param.Thrust * (t < t_burn);

    % Initialize derivative vector
    dY = zeros(4, 1);

    if t < param.t_turn
        dY(1) = T/m - D/m - g*sin(gamma0);
        dY(3) = -1/v * (g - v^2/(Re+h)) * cos(gamma0);
        dY(4) = Re/(Re+h) * v * cos(gamma0);
        dY(5) = v * sin(gamma0);
    else
        dY(1) = T/m - D/m - g*sin(gamma);
        dY(3) = -1/v * (g - v^2/(Re+h)) * cos(gamma);
        dY(4) = Re/(Re+h) * v * cos(gamma);
        dY(5) = v * sin(gamma);
    end

    % Prepare output struct for ode recall
    if nargout > 1
        parout.qdyn = qdyn;
        parout.Mach = norm([vx vy])/a;
    end

end