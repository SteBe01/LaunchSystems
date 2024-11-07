%% Ex1 - 25/10/2024

clear, clc
close all

% General parameters taken from lab
n = 7;                      % [-]         - Mass ratio m0/mf
I_sp = 390;                 % [s]         - Specific impulse
gamma0 = deg2rad(89.85);    % [rad]       - Initial pitch angle
param.m0 = 68000;           % [kg]        - Initial mass 
param.g0 = 9.81;            % [m/s^2]     - Earth's gravitational parameter
param.Re = 6378000;         % [m]         - Earth's radius
param.Thrust = 933.91e3;    % [N]         - Max Thrust
param.d = 5;                % [m]         - Rocket's diameter
param.TW = 1.4;             % [-]         - Thrust over Weight ratio
param.Cd = 0.5;             % [-]         - Drag coefficient
param.rho0 = 1.225;         % [kg/m^3]    - Air density at sea level
param.H0 = 7500;            % [m]         - Air density reference altitude

% Parameters taken with data from Pegasus
% param.m0 = 23130;           % [kg]        - Initial mass 
% param.d = 1.4;              % [m]         - Rocket's diameter


% Pitch maneuver
param.h_turn = 130;                       % [m]     - Initial maneuver altitude
param.gamma_turn = deg2rad(89.85);        % [rad]   - Initial flight path angle

% Compute necessary data
m_prop = param.m0*(1-1/n);        % [Kg]    - Propellant mass
param.m_dot = param.Thrust /(I_sp*param.g0);  % [Kg/s]  - Mass flow rate
param.t_burn = m_prop / param.m_dot;    % [s]     - Burning time

% ODE Initial state - V0 ~= 0 to avoid integration problems
y0 = [1 gamma0 0 0];
% options = odeset('RelTol',1e-8,'Events',@(x,y) event(x,y,1));
options = odeset('RelTol',1e-8);
[T,Y] = ode113(@(t, y) ballistic(t, y, param),[0 param.t_burn],y0,options);

% Retrieve data from ode
qdyn = zeros(1, length(T));
for ii = 1:length(T)
    [~, parout] = ballistic(T(ii), Y(ii, :), param);
    qdyn(ii) = parout.qdyn;
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



%% Functions

% function [value, isterminal, direction] = event(~,xx,isTerminal)
%     value = xx(4);
%     isterminal = isTerminal;
%     direction = 0;
% end

function [dY, parout] = ballistic(t,y, param)

    % Retrieve data from ode
    v = y(1);
    gamma = y(2);
    x = y(3);
    h = y(4);

    % Retrieve data used multiple times from param structure
    t_burn = param.t_burn;
    Re = param.Re;

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
    T = param.Thrust * (t < t_burn);

    % Initialize derivative vector
    dY = zeros(4, 1);

    if h < param.h_turn
        dY(1) = T/m - D/m - g;
        dY(4) = v;
    else
        dY(1) = T/m - D/m - g*sin(gamma);
        dY(2) = -1/v * (g - v^2/(Re+h)) * cos(gamma);
        dY(3) = Re/(Re+h) * v * cos(gamma);
        dY(4) = v * sin(gamma);
    end


    % if t < t_burn
    %     if h < param.h_turn
    %         dy = [T/m-D/m-g, 0, 0, v]';
    %     else
    %         dy = [T/m-D/m-g*sin(gamma)
    %             -1/v*(g-v^2/(Re+h))*cos(gamma)
    %             Re/(Re+h)*v*cos(gamma)
    %             v*sin(gamma)];
    %     end
    % else
    %     dy = [-D/m-g*sin(gamma)
    %         -1/v*(g-v^2/(Re+h))*cos(gamma)
    %         Re/(Re+h)*v*cos(gamma)
    %         v*sin(gamma)];
    % end

    % Prepare output struct for ode recall
    if nargout > 1
        parout.qdyn = qdyn;
    end

end