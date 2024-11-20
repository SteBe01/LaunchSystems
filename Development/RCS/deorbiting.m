% Parameters
clc
mu = 3.986e14; % [m^3/s^2] Gravitational parameter of Earth
Re = 6371e3;   % [m] Earth radius
rho0 = 1.225;  % [kg/m^3] Air density at sea level
H = 7e3;       % [m] Atmospheric scale height
Cd = 2.2;      % Drag coefficient
A = 1;         % [m^2] Satellite cross-sectional area
m = 500;       % [kg] Satellite mass
omega_earth = 7.2921159e-5; % [rad/s] Earth's angular velocity (rotational speed)

% Initial state (ECI coordinates)
r0 = [Re + 100e3; 0; 0]; % Initial position (600 km altitude)
v0 = [0; sqrt(mu / norm(r0)); 0]; % Circular velocity
state0 = [r0; v0]; % Initial state vector [r_x, r_y, r_z, v_x, v_y, v_z]

% Numerical integration
options = odeset('Events', @(t, state) impact_event(t, state, Re));    

tspan = [0, 1e5]; % Time span for simulation
[t, state] = ode45(@(t, state) orbital_dynamics(t, state, mu, Re, rho0, H, Cd, A, m, omega_earth), tspan, state0,options);

% Extract impact point
impact_position = state(end, 1:3); % Final position (ECI)
fprintf('Impact position (ECI): [%f, %f, %f] meters\n', impact_position);

% Convert to latitude and longitude
[lat, lon] = eci_to_geodetic(impact_position);
fprintf('Impact latitude: %.4f degrees\n', lat);
fprintf('Impact longitude: %.4f degrees\n', lon);
fprintf('Time to reach ground: %.4f hours\n', t(end)/3600);

% Plot trajectory
figure; 
plot3(state(:, 1) / 1e3, state(:, 2) / 1e3, state(:, 3) / 1e3, 'b');
hold on;
plot3(impact_position(1)/1000,impact_position(2)/1000,impact_position(3)/1000,'o')
plot3(0, 0, 0, 'ro', 'MarkerSize', 10); % Earth's center
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('Satellite Trajectory');
grid on;

axis equal;
hold on
x = @(x)cos(x)*Re/10^3; y = @(x)sin(x)*Re/10^3;
vec = linspace(0, 2*pi, 100);
plot3(x(vec), y(vec), zeros(100, 1)); % Earth's surface

% Orbital dynamics with drag and velocity relative to atmosphere
function dstate_dt = orbital_dynamics(t, state, mu, Re, rho0, H, Cd, A, m, omega_earth)
    % Unpack state vector
    r = state(1:3); % Position vector [x, y, z]
    v = state(4:6); % Velocity vector [vx, vy, vz]
    
    % Compute altitude
    h = norm(r) - Re; % Altitude above Earth's surface

    % Gravitational acceleration
    accel_gravity = -mu / norm(r)^3 * r;

    % Atmospheric density
    % if h > 0
    %     %rho = density_model(h);
    %     rho = rho0 * exp(-h / H); % Exponential atmospheric density decay with altitude
    % else
    %     rho = 0; % No atmosphere below surface
    % end
    % 
    rho = density_model(h);
    

    % Velocity of satellite relative to the atmosphere
    % Find the velocity of the atmosphere at the satellite's latitude
    lat = asin(r(3) / norm(r)); % Latitude from position vector
    v_atmosphere = omega_earth * Re * cos(lat); % Velocity of atmosphere at the satellite's latitude

    % Velocity of the satellite relative to the atmosphere
    v_relative = v - [v_atmosphere; 0; 0]; % Assuming atmosphere moves only in x direction (Eastward)

    % Drag acceleration (note: drag depends on relative velocity)
    accel_drag = -0.5 * Cd * A * rho * norm(v_relative) * v_relative / m;
       % Total acceleration
    accel = accel_gravity + accel_drag;

    % Derivatives of position and velocity
    dstate_dt = [v; accel];
end

% Event function to stop when satellite reaches Earth's surface
function [value, isterminal, direction] = impact_event(~, state, Re)
    % Stop condition: satellite reaches Earth's surface
    r = state(1:3); % Position vector [x, y, z]
    altitude = norm(r) - Re; % Altitude above Earth's surface
    
    % Trigger event when altitude < 10 meters
    value = altitude; % You can adjust this threshold
    isterminal = 1; % Stop integration
    direction = 0; % Only detect downward crossing
end

% Convert ECI position to geodetic coordinates (latitude, longitude)
function [lat, lon] = eci_to_geodetic(position)
    % Convert ECI position to geodetic coordinates (latitude, longitude)
    x = position(1);
    y = position(2);
    z = position(3);

    % Longitude
    lon = rad2deg(atan2(y, x));
    % Latitude
    lat = rad2deg(atan2(z, sqrt(x^2 + y^2)));
end
