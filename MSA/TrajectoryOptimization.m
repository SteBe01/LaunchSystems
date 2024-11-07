%% Ex1 - 25/10/2024

clear, clc
close all

g0 = 9.81;

m0 = 68000;                     % [kg]
n = 7;                          % m0/mf
I_sp = 390;                     % [s]
T = 933.91e3;                   % [N]
D = 5;                          % [m]
TW = 1.4;                       % T/W ratio
cd = 0.5;

% Solve for t_final or t_burn for the integration
% t_burn = m_prop / m_dot;
% m_dot = T /(I_sp*g0);
% m_prop = m0 - mf = m0*(1-1/n);

m_prop = m0*(1-1/n);
m_dot = T /(I_sp*g0);
t_burn = m_prop / m_dot;

h0 = 130;                       % [m]
gamma0 = deg2rad(89.85);        % [rad]

y0 = [1 gamma0 0 0];
t_vect = linspace(0, t_burn, 1000);
% options = odeset('RelTol',1e-8,'Events',@(x,y) event(x,y,1));
options = odeset('RelTol',1e-8);
[t,y] = ode113(@state_derivative,[0 t_burn],y0,options);

qdyn = zeros(1, length(t));

for ii = 1:length(t)
    [~, parout] = state_derivative(t(ii), y(ii, :));
    qdyn(ii) = parout.qdyn;
end


% Plots
boundary = 0;
subplot(2,2,1), hold on, grid on, title("v"), xlabel("Time [s]"), ylabel("Velocity [m/s]")
plot(t(1:end-boundary), y(1:end-boundary,1))
subplot(2,2,2), hold on, grid on, title("Gamma"), xlabel("Time [s]"), ylabel("Gamma [rad]")
plot(t(1:end-boundary), y(1:end-boundary,2))
subplot(2,2,3), hold on, grid on, title("Downrange"), xlabel("Time [s]"), ylabel("Downrange [m]")
plot(t(1:end-boundary), y(1:end-boundary,3))
subplot(2,2,4), hold on, grid on, title("h"), xlabel("Time [s]"), ylabel("Altitude [km]")
plot(t(1:end-boundary), y(1:end-boundary,4).*1e-3)

figure
plot(y(1:end, 4), qdyn./1e3);

figure
plot(y(1:end, 3)/1e3, y(1:end, 4)/1e3)



%% Functions

function [value, isterminal, direction] = event(~,xx,isTerminal)
    value = xx(4);
    isterminal = isTerminal;
    direction = 0;
end

function [dy, parout] = state_derivative(t,y)
Re = 6378000;
g0 = 9.81;

m0 = 68000;                     % [kg]
n = 7;                          % m0/mf
I_sp = 390;                     % [s]
T = 933.91e3;                   % [N]
D = 5;                          % [m]
TW = 1.4;                       % T/W ratio
cd = 0.5;

h0 = 130;                       % [m]

v = y(1);
gamma = y(2);
x = y(3);
h = y(4);

S = pi*(D^2/4);

hturn = 130;
rho0 = 1.225;
h0 = 7500;
rho = @(hv) rho0*exp(-h/h0);
m_prop = m0*(1-1/n);
m_dot = T /(I_sp*g0);
t_burn = m_prop / m_dot;

m = @(t) m0 - ((m_dot*t)*(t<t_burn) + (m_dot*t_burn)*(t>=t_burn));

g = @(xx) g0/((1+xx/Re)^2);

qdyn = 0.5*rho(h)*v^2;
D = qdyn*S*cd;


if t < t_burn
    if h < hturn
        dy = [T/m(t)-D/m(t)-g(h), 0, 0, v]';
    else
        dy = [T/m(t)-D/m(t)-g(h)*sin(gamma)
        -1/v*(g(h)-v^2/(Re+h))*cos(gamma)
        Re/(Re+h)*v*cos(gamma)
        v*sin(gamma)];
    end
else
    dy = [-D/m(t)-g(h)*sin(gamma)
    -1/v*(g(h)-v^2/(Re+h))*cos(gamma)
    Re/(Re+h)*v*cos(gamma)
    v*sin(gamma)];
end

if nargout > 1
    parout.qdyn = qdyn;
end

end