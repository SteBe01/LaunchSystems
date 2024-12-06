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
y0_stg1 = [0 0 0, 0 20 -100, 0.01 0.01 0.01, 1 0 0 0]';

options_stg1 = odeset('RelTol',1e-8, 'MaxStep', 0.1, 'Events', @(t, y) touchdown(t, y));
[tvect, yvect] = ode113(@(t,y) dyn(t, y, stages.stg1, params), [0 1e3], y0_stg1, options_stg1);

T = tvect;
Y = yvect;

downrange = params.Re./(params.Re+Y(:,1)) .* Y(:,1);

% dcm = zeros(3,3,length(T));
% for ii = 1:length(T)
%     [~, parout] = dyn(T(ii),Y(ii,:),stages.stg1,params);
%     dcm(:,:,ii) = parout.dcm';
% end

%% Animation

% animateOrientation(dcm, T, "Palle1");



%% Plots for simple 3 dof

set(0, 'DefaultLineLineWidth', 1.5)

figure, hold on, axis equal, grid on
plot3(Y(:,1),Y(:,2),-Y(:,3))
view(40,30)

%%

% figure, hold on
% plot(T, Y(:,7))
% plot(T, Y(:,8))
% plot(T, Y(:,9))




%% Functions

function [value, isterminal, direction] = touchdown(~, y)
    value = y(3);
    isterminal = 1;
    direction = 0;
end

