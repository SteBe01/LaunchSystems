%% Lab 18/11/2024

clear, clc
close all

p.Ma = 3.75;                % 1/s1
p.Md = 4.54;                % 1/s2
p.F = 1.6681e6;             % [N]
p.m = 85082.45;             % [kg]
p.V = 402.336;              % [m/s]
p.N = 1.0676e6;
p.T = 1.5168e6;             % [N]
p.alphawdot = 0;

alphaw = 5.73*pi/180;       % [rad]

u = 0;

x0 = [0 0 alphaw 0 0]';

[tvect, xvect] = ode45(@(t,x) fun_dyn(t,x,p,u), [0 5], x0);

subplot(4,1,1), hold on, title("No control")
plot(tvect, xvect(:,1).*180/pi, 'k', 'LineWidth',2)
grid on, ylabel('\Theta [deg]')
subplot(4,1,2)
plot(tvect, u*ones(length(tvect)), 'k', 'LineWidth',2)
grid on, ylabel('\delta [deg]')
subplot(4,1,3)
plot(tvect, xvect(:,3).*180/pi, 'k', 'LineWidth',2)
grid on, ylabel('\alpha [deg]')
subplot(4,1,4)
plot(tvect, xvect(:,4), 'k', 'LineWidth',2)
grid on, ylabel('z [m]')
xlabel('Time [s]')


%% PD control
% delta = u = -K1*THETA - K2*DOTTHETA

k1 = 2;
k2 = 0.8;

h = 0.01;
tspan = 0:h:5;

x0 = [0 0 alphaw 0 0]';
X = zeros(length(tspan),length(x0));
X(1,:) = x0';

u_c = zeros(length(tspan)-1,1);
for ii = 1:length(tspan)-1
    if ii == 1
        u_c(ii) = 0;
    else
        u_c(ii) = -k1*x0(1) - k2*x0(2);
    end
    [tvect_c, xvect_c] = ode45(@(t,x) fun_dyn(t,x,p,u_c(ii)), [tspan(ii) tspan(ii+1)], x0);
    x0 = xvect_c(end,:)';
    X(ii+1,:) = xvect_c(end,:);
end

figure
subplot(4,1,1), hold on, title("Control")
plot(tspan, X(:,1).*180/pi, 'k', 'LineWidth',2)
grid on, ylabel('\Theta [deg]')
subplot(4,1,2)
plot(tspan(2:end), u_c.*180/pi, 'k', 'LineWidth',2)
grid on, ylabel('\delta [deg]')
subplot(4,1,3)
plot(tspan, X(:,3).*180/pi, 'k', 'LineWidth',2)
grid on, ylabel('\alpha [deg]')
subplot(4,1,4)
plot(tspan, X(:,4), 'k', 'LineWidth',2)
grid on, ylabel('z [m]')
xlabel('Time [s]')


%% Feedback Drift Minimum

k1 = 2;
k2 = 0.8;
k3 = 3.614;

h = 0.01;
tspan = 0:h:5;

x0 = [0 0 alphaw 0 0]';
X = zeros(length(tspan),length(x0));
X(1,:) = x0';

u_c = zeros(length(tspan)-1,1);
for ii = 1:length(tspan)-1
    if ii == 1
        u_c(ii) = 0;
    else
        u_c(ii) = -k1*x0(1) - k2*x0(2) - k3*x0(3);
    end
    [tvect_c, xvect_c] = ode45(@(t,x) fun_dyn(t,x,p,u_c(ii)), [tspan(ii) tspan(ii+1)], x0);
    x0 = xvect_c(end,:)';
    X(ii+1,:) = xvect_c(end,:);
end

figure
subplot(4,1,1), hold on, title("Feedback Drift Minimum")
plot(tspan, X(:,1).*180/pi, 'k', 'LineWidth',2)
grid on, ylabel('\Theta [deg]')
subplot(4,1,2)
plot(tspan(2:end), u_c.*180/pi, 'k', 'LineWidth',2)
grid on, ylabel('\delta [deg]')
subplot(4,1,3)
plot(tspan, X(:,3).*180/pi, 'k', 'LineWidth',2)
grid on, ylabel('\alpha [deg]')
subplot(4,1,4)
plot(tspan, X(:,4), 'k', 'LineWidth',2)
grid on, ylabel('z [m]')
xlabel('Time [s]')


%% Feedback Load Minimum

k1 = 0;
k2 = 0.8;
k3 = 3.614;

h = 0.01;
tspan = 0:h:5;

x0 = [0 0 alphaw 0 0]';
X = zeros(length(tspan),length(x0));
X(1,:) = x0';

u_c = zeros(length(tspan)-1,1);
for ii = 1:length(tspan)-1
    if ii == 1
        u_c(ii) = 0;
    else
        u_c(ii) = -k1*x0(1) - k2*x0(2) - k3*x0(3);
    end
    [tvect_c, xvect_c] = ode45(@(t,x) fun_dyn(t,x,p,u_c(ii)), [tspan(ii) tspan(ii+1)], x0);
    x0 = xvect_c(end,:)';
    X(ii+1,:) = xvect_c(end,:);
end

figure
subplot(4,1,1), hold on, title("Feedback Load Minimum")
plot(tspan, X(:,1).*180/pi, 'k', 'LineWidth',2)
grid on, ylabel('\Theta [deg]')
subplot(4,1,2)
plot(tspan(2:end), u_c.*180/pi, 'k', 'LineWidth',2)
grid on, ylabel('\delta [deg]')
subplot(4,1,3)
plot(tspan, X(:,3).*180/pi, 'k', 'LineWidth',2)
grid on, ylabel('\alpha [deg]')
subplot(4,1,4)
plot(tspan, X(:,4), 'k', 'LineWidth',2)
grid on, ylabel('z [m]')
xlabel('Time [s]')

