function [LOAD,TR] = Attachments(TR,loads)



g0 = 9.80665; %m/s^2
FS=TR.FS;
mat_id_1 = 2;
m_tot = TR.m_tot;

[MAT_1] = material_selection(mat_id_1);
sigma_yield = MAT_1.sy;

x_com=TR.x_min;
nx=loads.nx;
nz=loads.nz;
a=TR.a;
b=TR.b;
l=TR.l;
J=TR.J;
R_1 = TR.Diam1/2;
F_z = (nz+1)*g0*m_tot;
F_x = (nx)*g0*m_tot;
LOAD.F_z = F_z;
R_a = F_z*0.5;
R_b = F_z*0.5;
M=R_a*x_com;
P_a = F_x;
thick=TR.t;
sigma_a = FS*abs(P_a/(2*pi*(R_1)));
sigma_b = FS*(M/(pi*R_1^2));
sigma_max = sigma_a + sigma_b;
LOAD.sigma_max = sigma_max;
LOAD.thick_2 = sigma_max/sigma_yield;

if LOAD.thick_2<thick

    fprintf('Correct dimensioning');

else 
 fprintf('Wrong dimensioning!!!!!');


end

%% Plot:

% Segment lengths
x_segments = [x_com,l];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting

% Axial Forces (P) at segment ends
P_forces = [0,0];
LOAD.P_forces = [P_a,0,0];

figure;
hold on;
% for i = 1:length(x_segments)
%     plot([x(i+1), x(i+1)], [P_a, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
% end
stairs(x, [P_a, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Axial Force [N]');
title('Axial Force Diagram for Worst Case Carrier Maneuver');
xlim([0, l]);
grid on;

% Segment lengths
x_segments = [x_com,l];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting

% Axial Forces (P) at segment ends
P_forces = -[R_b,R_b];
figure;
hold on;
% for i = 1:length(x_segments)
%      plot([x(i+1), x(i+1)], [R_a, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
%  end
stairs(x, [R_a, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Shear Force [N]');
title('Shear Force Diagram for Worst Case Carrier Maneuver');
xlim([0, l]);
grid on;


% Plot the bending moment diagram
figure;
hold on;

plot([0,a,l],[0, -M,0], 'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Bending Moment [Nm]');
title('Bending Moment Diagram for Worst Case Carrier Maneuver');
xlim([0, l]);
grid on;


end