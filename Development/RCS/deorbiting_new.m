% Parameters
clc
par.mu = 398600; % [km^3/s^2] Gravitational parameter of Earth
par.Re = 6371;   % [km] Earth radius
par.rho0 = 1.225;  % [kg/m^3] Air density at sea level
par.H = 7e3;       % [m] Atmospheric scale height
geom.Cd = 2.2;      % Drag coefficient
geom.A = 0.95;         % [m^2] cross-sectional area
geom.m = 150;       % [kg] mass
par.omega_earth = 7.2921159e-5; % [rad/s] Earth's angular velocity (rotational speed)
h0 = 400;  %[km]

% Initial state (ECI coordinates)
r0 = [par.Re + h0; 0; 0];
v0 = [0; sqrt(par.mu / norm(r0)); 0]; % Circular velocity
[r0,v0] = kep2car([norm(r0),0,deg2rad(90),0,0,deg2rad(30)],par.mu);
state0 = [r0; v0]; % Initial state vector [r_x, r_y, r_z, v_x, v_y, v_z]


tspan = [0, 10^8]; % Time span for simulation
options = odeset('Events', @(t,y) impact_event(t, y, par.Re));   

% [~, y] = ode78(@(t, y)ode_2bp(t,y,par.mu,par,geom), tspan, state0,options);
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 ,'Events', @(t,y) impact_event(t, y, par.Re));
[ t, Y ] = ode113( @(t,y) ode_2bp_perturbed( t, y, par.mu, par.Re, 0, par.omega_earth, geom.A, geom.Cd), tspan, state0, options );
figure
plot3(Y(:,1),Y(:,2),Y(:,3),LineWidth=1.5)
axis equal;
hold on
Earth_3D(par.Re)
plot3(Y(end,1),Y(end,2),Y(end,3),'o')


% Extract impact point
impact_position = Y(end, 1:3); % Final position (ECI)
fprintf('Impact position (ECI): [%f, %f, %f] meters\n', impact_position);

% Convert to latitude and longitude
[lat, lon] = eci_to_geodetic(impact_position);
fprintf('Impact latitude: %.4f degrees\n', lat);
fprintf('Impact longitude: %.4f degrees\n', lon);
fprintf('Time to reach ground: %.4f hours\n', t(end)/3600);

function dy = ode_2bp(~, y,mu, par,geom)

r = y(1:3);
v = y(4:6);
mu = par.mu;
h = norm(r) - par.Re; % Altitude above Earth's surface

rho = density_model(h);

v_atmosphere = cross([0 0 par.omega_earth]',r);

v_relative = v - v_atmosphere; % Assuming atmosphere moves only in x direction (Eastward)

accel_drag = -0.5 * geom.Cd * geom.A * rho * norm(v_relative) * v_relative / geom.m;

dy = [v; (-mu/(norm(r)^3)) * r + accel_drag];

end

% Event function to stop when satellite reaches Earth's surface
function [value, isterminal, direction] = impact_event(~, state, Re)
    r = state(1:3); % Position vector [x, y, z]
    altitude = norm(r) - Re; % Altitude above Earth's surface
    value = altitude-0.001; 
    if altitude <=0
        value = 0;
    end
    
    isterminal = 1; % Stop integration
    direction = 0; % versus (0 detects both)
end

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
