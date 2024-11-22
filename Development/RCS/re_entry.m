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
geom.R_n = 0.5;           % Nose radius of the launcher (km)
par.emissivity = 0.9;        % Emissivity of the surface
par.sigmaSB = 5.67e-8;      % Stefan-Boltzmann constant (W/m^2/K^4, converted later)
par.theta_g_0 = deg2rad(0);          % Initial Greenwich meridian position         [rad]

% Initial state (ECI coordinates)
r0 = [par.Re + h0; 0; 0];
[r0,v0] = kep2car([norm(r0),0,deg2rad(98),0,0,deg2rad(180)],par.mu);
state0 = [r0; v0; T0];


tspan = [0, 10^8]; % Time span for simulation

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 ,'Events', @(t,y) impact_event(t, y, par.Re));
[ t, Y ] = ode113( @(t,y) ode_2bp_drag( t, y,par,geom), tspan, state0, options );

% Extract impact point
impact_position = Y(end, 1:3); % Final position (ECI)
fprintf('Impact position (ECI): [%f, %f, %f] meters\n', impact_position);

% Convert to latitude and longitude
[lat, lon] = eci_to_geodetic(impact_position);
fprintf('Impact latitude: %.4f degrees\n', lat);
fprintf('Impact longitude: %.4f degrees\n', lon);
fprintf('Time to reach ground: %.4f hours\n', t(end)/3600);

radius = zeros(length(Y(:,1)),1);
for i = 1:length(Y(:,1))
    radius(i) = norm(Y(i,1:3)) -par.Re;  
end

q_net = heat_flux(Y,geom,par);

%%%% Plots
figure
plot3(Y(:,1),Y(:,2),Y(:,3),LineWidth=1.5)
axis equal;
hold on
Earth_3D(par.Re)
plot3(Y(end,1),Y(end,2),Y(end,3),'o')

figure
plot(radius,Y(:,4:6))
legend('V_x','V_y','V_z'); xlabel('Altitude [km]'); ylabel('Velocity [km/s]');grid on;set ( gca, 'XDir', 'reverse' ) 

figure
plot(radius,Y(:,7))
xlabel('Altitude [km]'); ylabel('Temperature peak [K]');grid on
set ( gca, 'XDir', 'reverse' )

figure
plot(radius,q_net)
xlabel('Altitude [km]'); ylabel('Heat flux [W/m^2]');grid on
set ( gca, 'XDir', 'reverse' )

ground_track(t,Y,par)

%%% Functions
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

    lon = rad2deg(atan2(y, x));               % Longitude
    lat = rad2deg(atan2(z, sqrt(x^2 + y^2))); % Latitude
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

a_drag = (-0.5 * geom.A/geom.m * rho * geom.Cd * norm(v_rel)^2) * (v_rel ./ norm(v_rel)) * 1e-3;        % [m/s^2]

% Set the derivatives of the state
dy = [ v;
           (-par.mu/(norm(r))^3)*r + a_drag];

% Temperature
Q_aero = par.kappa * sqrt(rho/geom.R_n) * (norm(v_rel))^3*geom.A;  % Convert v to m/s for heat flux POTREBBE MANCARE L'AREA
% q_solar = (1-par.alpha_alb)*par.xhi*par.q_sol*geom.A; 
Q_space = par.emissivity*par.sigmaSB*y(7)^4*geom.A;
Q_net = Q_aero -Q_space;
dy(7) = Q_net / (geom.m * par.cp);       % Temperature rate (K/s)
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

function ground_track(t_vec,y_vec,par)
delta = zeros(length(t_vec), 1);
alpha = zeros(length(t_vec), 1);
lat   = zeros(length(t_vec), 1);
lon   = zeros(length(t_vec), 1);

for i = 1:length(t_vec)

    x = real(y_vec(i, 1));
    y = real(y_vec(i, 2));
    z = real(y_vec(i, 3));

    r = norm([x, y, z]);
    
    delta(i) = asin(z/r);
    alpha(i) = atan2(y, x);

    % Conversion to longitude and latitude
    theta_g = par.theta_g_0 + par.om_E*(t_vec(i) - t_vec(1));

    lon(i) = rad2deg(wrapToPi(alpha(i) - theta_g));
    lat(i) = rad2deg((delta(i)));
end
   
    figure('Position', [10 10 720 360])
    plot(lon, lat, 'r', 'LineWidth', 1, 'HandleVisibility','off');
    hold on

    axis([-180 180 -90 90])
    h = image(xlim, -ylim, imread('earth_texture_2D.jpg')); 
    uistack(h, 'bottom')
    
    % initial and final points
    plot(lon(1), lat(1), 'g*', 'LineWidth', 1, 'MarkerSize', 5, 'DisplayName', 'Start')
    hold on
    plot(lon(end), lat(end), 'g^', 'LineWidth', 1, 'MarkerSize', 5, 'DisplayName', 'End')
    hold on
    
    xlabel('Longitude [deg]'), ylabel('Latitude [deg]')
    legend show
   
end

function q_net = heat_flux(Y,geom,par)
q_net = zeros(length(Y(:,1)),1);
for i = 1:length(Y(:,1))
r = Y(i,1:3);
v = Y(i,4:6);
h = norm(r)- par.Re;
rho = density_model(h);
v_rel = v - cross([0 0 par.om_E], r)';          % [km]
v_rel = v_rel * 1e3;                        % [m/s]
q_aero = par.kappa * sqrt(rho/geom.R_n) * (norm(v_rel))^3;  % Convert v to m/s for heat flux POTREBBE MANCARE L'AREA
q_space = par.emissivity*par.sigmaSB*Y(i,7)^4;
q_net(i) = q_aero -q_space;
end
end