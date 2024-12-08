function [CLAMP] = Attachments_FINAL(GEOMETRY,LOAD,PlotFlag,SaveFlag)


g0 = 9.80665; %m/s^2
m2h =GEOMETRY.m1+GEOMETRY.m2+GEOMETRY.m3+GEOMETRY.m4+GEOMETRY.m5;
m1h =GEOMETRY.m6+GEOMETRY.m7+GEOMETRY.m8;
m3h = GEOMETRY.m9+GEOMETRY.m10;
Diam_1 = GEOMETRY.Diam_1;
Diam_2 = GEOMETRY.Diam_2;
t=max([GEOMETRY.thick_1_max;GEOMETRY.thick_2_max]);
Rad_2 = Diam_1/2;
Rad_1 = Rad_2 - t;
L_tot = GEOMETRY.l_tot;
Cd_c = LOAD.Cd;
Cl_c = LOAD.Cl;
q = 0.5*1.1*(0.6*343)^2;
%q = LOAD.q;
D = q*(pi*Rad_2^2)*Cd_c;
nx_c = LOAD.nx_c;
nz_c = LOAD.nz_c;
F1_z = m2h*g0*(nz_c+1);
F2_z = m1h*g0*(nz_c+1);
F1_x = m2h*g0*(nx_c);
F2_x = m1h*g0*(nx_c);
F3_x = m3h*g0*(nx_c);
F3_z = m2h*g0*(nz_c+1);
[X_COM2] =MASS_PROPERTIES_2(GEOMETRY.m_prop_2,GEOMETRY);

syms L_cone F1 F2 F3 L_fin R1 R2 a b c d e z L

eq1 = L_cone -F1 +R1-F2 -F3+L_fin+R2==0;
eq2 = -L_cone*z + F1*a -F2*b + R2*(b+c) -F3*(d+c+b) +L_fin*(e+d+c+b)==0;
eq3 = -L_cone*(z+a+b+c)+F1*(a+b+c)-R1*(b+c)+F2*c -F3*d+L_fin*(e+d)==0;
eq4 = L_cone+L_fin-L==0;

sol = solve([eq1, eq2, eq3,eq4], [R1, R2,L_cone,L_fin]);

F1=F1_z;
F2=F2_z;
F3=F3_z;
L=q*(pi*Rad_2^2)*Cl_c;
param=0.467;
L_cone = L*param;
L_fin = L*(1-param);
l_cone = GEOMETRY.b1*2/3;
z = X_COM2-l_cone;
a =  GEOMETRY.b1+GEOMETRY.b2+GEOMETRY.b3+GEOMETRY.b34+GEOMETRY.b4+GEOMETRY.b5+GEOMETRY.b6 - X_COM2;
l_constraint = 8;
b = 2*l_constraint/3;
c = l_constraint/3;
d = (L_tot-(l_cone+z+a+b+c))/2;
X_COM3 = l_cone+z+a+b+c+d;
e =d;
CHECK=l_cone+z+a+b+c+d+e-L_tot;
if CHECK>10^-3

fprintf('ERROR IN LENGTH');

end

P1 =  D + F1_x+F2_x+F3_x;
CLAMP.P1=P1;

sol = subs(sol);


if PlotFlag==1

x_segments = [l_cone+z,a,b,c,d,e];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting

% Axial Forces (P) at segment ends
P_forces = [F1_x+D,F1_x+D-P1,F1_x+D-P1+F2_x,F1_x+D-P1+F2_x,0,0];

figure()
tiledlayout(3,1);
nexttile;
hold on;
for i = 1:length(x_segments)
     plot([x(i+1), x(i+1)], [D, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
end
stairs(x, [D, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Axial Force [N]');
title('Axial Force Diagram for Worst Case Carrier Maneuver');
xlim([0,L_tot]);
grid on;


R1 = double(subs(sol.R1));
R2 = double(subs(sol.R2));

CLAMP.R1=R1;
CLAMP.R2=R2;
x_segments = [l_cone,z,a,b,c,d,e];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting

% Axial Forces (P) at segment ends
%P_forces = [L_cone,L_cone-F1,L_cone-F1+R1,L_cone-F1+R1-F2,L_cone-F1+R1-F2+R2,L_cone-F1+R1-F2+R2-F3,0];
P_forces = [-(L_cone-(F1*2/3)),L_cone-(F1*2/3)-F1,L_cone-F1+R1,L_cone-F1+R1-F2,L_cone-F1+R1-F2+R2,L_cone-F1+R1-F2+R2-F3,0];

nexttile; % SHEAR
hold on;
for i = 1:length(x_segments)
    plot([x(i+1), x(i+1)], [0, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
end
stairs(x, [0, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Shear Force [N]');
title('Shear Force Diagram for Worst Case Carrier Maneuver');
xlim([0,L_tot]);
grid on;



P_forces = [L_cone,L_cone-F1,-(L_cone-F1+R1),-(L_cone-F1+R1-F2),-(L_cone-F1+R1-F2+R2),(L_cone-F1+R1-F2+R2-F3),0];
x_segments = [l_cone,z,a,b,c,d,e];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting
    % Calculate the bending moment
M_moments = zeros(size(P_forces));
M_moments(1) = 0; % Initial bending moment

for i = 2:length(P_forces)
    M_moments(i) = M_moments(i-1) + P_forces(i-1) * x_segments(i-1);
end

M_moments(2)=-L_cone*z;
M_moments(3)=0;

M_moments(4) = min(M_moments)*0.7;

M_moments(5)=0;

M_moments(end-1)=-M_moments(end-1);

M_moments(end)=0;


% Plot the bending moment diagram
nexttile;
hold on;
for i = 1:length(x_segments)
    plot([x(i+1), x(i+1)], [0, M_moments(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
end
plot(x,[0, M_moments], 'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Bending Moment [Nm]');
title('Bending Moment Diagram for Worst Case Carrier Maneuver');
xlim([0, sum(x_segments)]);
grid on;
%yticks = get(gca, 'YTick');                         
% yticklabels = arrayfun(@(v) sprintf('%d \\times 10^5', abs(v) / 1e5), yticks, 'UniformOutput', false);
% set(gca, 'YTickLabel', yticklabels);
yticks = get(gca, 'YTick'); % Get current tick values

% Convert tick values to scientific notation strings with "e" notation and absolute values
yticklabels = arrayfun(@(y) sprintf('%0.2e', abs(y)), yticks, 'UniformOutput', false);

% Set the formatted tick labels to the y-axis
set(gca, 'YTickLabel', yticklabels);

end

if SaveFlag==1
exportgraphics(gcf, 'Load Diagram Attachments.pdf', 'ContentType','vector');
end

end % function