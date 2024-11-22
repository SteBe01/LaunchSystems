% Parameters
clc
par.mu = 398600; % [km^3/s^2] Gravitational parameter of Earth
par.Re = 6371;   % [km] Earth radius
geom.Cd = 2.2;      % Drag coefficient
geom.A = 0.95;         % [m^2] cross-sectional area
geom.m = 150;       % [kg] mass
par.om_E = 7.2921159e-5; % [rad/s] Earth's angular velocity (rotational speed)
h0 = 100;  %[km]
T0 = 300; %[k]
par.kappa = 1.83e-4;  % Heat transfer constant (W/km^2/K^0.5)
par.cp = 900;  % Specific heat capacity (J/kg/K)
par.thermal_mass = 150;  % Thermal mass of the vehicle in thermal equilibrium (kg)


% Initial state (ECI coordinates)
r0 = [par.Re + h0; 0; 0];
v0 = [0; sqrt(par.mu / norm(r0)); 0]; % Circular velocity
[r0,v0] = kep2car([norm(r0),0,deg2rad(98),0,0,deg2rad(30)],par.mu);
state0 = [r0; v0;T0]; % Initial state vector [r_x, r_y, r_z, v_x, v_y, v_z]


tspan = [0, 10^8]; % Time span for simulation
options = odeset('Events', @(t,y) impact_event(t, y, par.Re));   

% [~, y] = ode78(@(t, y)ode_2bp(t,y,par.mu,par,geom), tspan, state0,options);
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 ,'Events', @(t,y) impact_event(t, y, par.Re));
[ t, Y ] = ode113( @(t,y) ode_2bp_drag( t, y,par,geom), tspan, state0, options );
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

% Event function to stop when satellite reaches Earth's surface
function [value, isterminal, direction] = impact_event(~, state, Re)
    r = state(1:3); % Position vector [x, y, z]
    altitude = norm(r) - Re; % Altitude above Earth's surface
    value = altitude-0.1; 
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

function dy = ode_2bp_drag( ~, y,par,geom)

r = y(1:3);
v = y(4:6);


h = norm(r) - par.Re ;                        
if h <= 0
    error("Impact on Earth detected! Reduce simulation time");
end

rho = density_model(h);                         %[kg/m^3]
 
v_rel = v - cross([0 0 par.om_E], r)';          % [km]
v_rel = v_rel * 1e3;                        % [m/s]

a_drag = (-0.5 * geom.A/geom.m * rho * geom.Cd * norm(v_rel)^2) ...
    * (v_rel ./ norm(v_rel)) * 1e-3;        % [m/s^2]

% Set the derivatives of the state
dy = [ v
 (-par.mu/(norm(r))^3)*r + a_drag];

% Temperature
 q_aero = par.kappa * sqrt(rho) * (norm(v) * 1e3)^3;  % Convert v to m/s for heat flux POTREBBE MANCARE L'AREA
% q_solar = (1-par.alpha_alb)*par.xhi*par.q_sol*geom.A;  %ricordati A e non A/m
% q_space = par.emissivity*par.sigmaSB*y(7)^4*par.A;
 q_net = q_aero;
dy(7) = q_net / (par.thermal_mass * par.cp);       % Temperature rate (K/s)


end

function rho = density_model(h)

%input  h: altitude in [km]

h0_vect = [0 25 30 40 50 60 70 80 90 100 110 120 130 140 150 ...
        180 200 250 300 350 400 450 500 600 700 800 900 1000]';
    rho0_vect = [1.225 3.899*1e-2 1.774*1e-2 3.972*1e-3 1.057*1e-3 ...
        3.206*1e-4 8.770*1e-5 1.905*1e-5 3.396*1e-6 5.297*1e-7 9.661*1e-8 ...
        2.438*1e-8 8.484*1e-9 3.845*1e-9 2.070*1e-9 5.464*1e-10 ...
        2.789*1e-10 7.248*1e-11 2.418*1e-11 9.158*1e-12 3.725*1e-12 ...
        1.585*1e-12 6.967*1e-13 1.454*1e-13 3.614*1e-14 1.170*1e-14 ...
        5.245*1e-15 3.019*1e-15]';
    H_vect = [7.249 6.349 6.682 7.554 8.382 7.714 6.549 5.799 5.382 5.877 ...
        7.263 9.473 12.636 16.149 22.523 29.740 37.105 45.546 53.628 ...
        53.298 58.515 60.828 63.822 71.835 88.667 124.64 181.05 268.00]';

if h < 0
    rho0 =0;
    h0 = 0;                            % [km]
    H = 0;                              % [km]
else
    idx = sum(h0_vect < h);

    rho0 = rho0_vect(idx);
    h0 = h0_vect(idx);                            % [km]
    H = H_vect(idx);                              % [km]
end
    rho = rho0 * exp(-(h-h0)/H);                % [kg/m^3]

end


