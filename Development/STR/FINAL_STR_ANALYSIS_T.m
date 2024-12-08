function [STRUCT]=FINAL_STR_ANALYSIS_T(GEOMETRY,FORCES,CASE,SaveFlag)

b1 = GEOMETRY.b1; % fairing length
b2 = GEOMETRY.b2; % forward skirt length
b3 = GEOMETRY.b3; % height of fuel tank2
b34 = GEOMETRY.b34; % 16 cm between tanks
b4 = GEOMETRY.b4; %  height of oxygen tank2
b5 = GEOMETRY.b5; %  aft skirt length 2
b6 = GEOMETRY.b6; % interstage length
b7 = GEOMETRY.b7; % height of fuel tank1
b78 = GEOMETRY.b78; % 16 cm between tanks
b8 = GEOMETRY.b8; % height of oxygen tank1
b9 = GEOMETRY.b9; % aft skirt length 1
b10 = GEOMETRY.b10; % h engine

x_segments = [b1, b2, b3,b34 ,b4,b5,b6,b7,b78,b8,b9,b10];
l_tot=sum(x_segments);

m1 = GEOMETRY.m1; 
m2 = GEOMETRY.m2;
m3 = GEOMETRY.m3;
m4 = GEOMETRY.m4;
m5 = GEOMETRY.m5; 
m6 = GEOMETRY.m6;
m7 = GEOMETRY.m7;
m8 = GEOMETRY.m8;
m9 = GEOMETRY.m9;
m10 = GEOMETRY.m10; 

[x_com] = MASS_PROPERTIES_1(GEOMETRY.m_prop_1,GEOMETRY);
[x_com2] = MASS_PROPERTIES_2(GEOMETRY.m_prop_2,GEOMETRY);
m=m1+m2+m3+m4+m5+m6+m7+m8+m9+m10;

g0 = 9.80665; %m/s^2

q = FORCES.q;
A_ref = pi*(GEOMETRY.Diam_1^2)/4;
alpha = FORCES.alpha;

AER=Aer_force(q,alpha);

Cl_nose = FORCES.Cl;
Cl_fin = FORCES.Cl;
Cd = FORCES.Cd;

D = q*Cd*A_ref;
L_fin = q*Cl_fin*A_ref;
L_nose = q*Cl_nose*A_ref;

delta = FORCES.delta;


switch CASE % 1: 1 piece; 2: 10 pieces; 3: 4 pieces

    case 1
STRUCT=1;

%T = FORCES.T;

%nx =((T*cos(delta) - D)/(m*g0));
%nx =((T*cos(delta) - D - m*g0)/(m*g0));

 nx = FORCES.nx;
% 
% D = -nx*m*g0 +T*cos(delta);
T = (D + nx*m*g0)/(cos(delta));
nz=FORCES.nz;

L = -T*sin(delta) + m*g0 + nz*g0*m;
L_cone = (0.45)*L;
L_fin = (0.55)*L;
arm = (nz*g0*m*(l_tot-x_com))/L_cone;
%arm = ((nz+1)*g0*m*(l_tot-x_com))/L_cone;

mra = l_tot-arm;
%nz = FORCES.nz;

%% Load:

%% 1.
D1 = D;

P1 = D1;
R1 = -L_cone;

%% 2.
P2 = P1 + m*g0*(nx);
R2 = R1 + m*g0*(nz+1);
%% 3. 

P3 = P2  - T*cos(delta);
R3 = R2 - T*sin(delta)-L_fin;

% Moments:

%M_1 = 

M_cg = -L_cone*(x_com-mra) + (L_fin+T*sin(delta))*(l_tot-x_com);

%% PLOTS:

% Segment lengths
x_segments = [x_com,l_tot-x_com];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting

% Axial Forces (P) at segment ends
P_forces = [P2,P2];

figure();
tiledlayout(3,1)
nexttile
hold on;
for i = 1:length(x_segments)
    plot([x(i+1), x(i+1)], [D, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
end
stairs(x, [D, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Axial Force [N]');
title('Axial Force Diagram for Max-Q');
xlim([0, l_tot]);
grid on;

% Segment lengths
x_segments = [b1*2/3,x_com-(b1*2/3),l_tot-x_com];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting

% Axial Forces (P) at segment ends
P_forces = -[R1,R2,R2];
nexttile
hold on;
for i = 1:length(x_segments)
     plot([x(i+1), x(i+1)], [0, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
 end
stairs(x, [0, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Shear Force [N]');
title('Shear Force Diagram for Max-Q');
xlim([0, sum(x_segments)]);
grid on;

P_forces = [R1,R2,R2];
x_segments =[b1*2/3,x_com-(b1*2/3),l_tot-x_com] ;
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting
    % Calculate the bending moment
M_moments = zeros(size(P_forces));
M_moments(1) = 0; % Initial bending moment

for i = 2:length(P_forces)
    M_moments(i) = M_moments(i-1) + P_forces(i-1) * x_segments(i-1);
end

M_moments(end)=0;

x_inter = [b1*2/3;x_com];
y_inter = [0;M_moments(2)];
M_inter_eq = polyfit(x_inter,y_inter,1);
M_inter_line = @(x) M_inter_eq(1)*x + M_inter_eq(2);

x_M_2 = [0,0,sum(b1+b2+b3+b4+b5+b34),sum(b1+b2+b3+b4+b5+b34),0];
y_M_2 = [0,M_inter_line(sum(b1+b2+b3+b4+b5+b34)),M_inter_line(sum(b1+b2+b3+b4+b5+b34)),0,0];
x_M_1 = [sum(b1+b2+b3+b4+b5+b34),sum(b1+b2+b3+b4+b5+b34),l_tot,l_tot,sum(b1+b2+b3+b4+b5)];
y_M_1 = [0,min(M_moments),min(M_moments),0,0];

% Plot the bending moment diagram
nexttile
hold on;
for i = 1:length(x_segments)
    plot([x(i+1), x(i+1)], [0, M_moments(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
end
plot(x,[0, M_moments], 'LineWidth', 2, 'Color', 'b'); % Step-style plot
hold on;
plot(x_M_2,y_M_2,'Color','r','LineStyle','--','LineWidth', 1);
hold on;
plot(x_M_1,y_M_1,'Color','r','LineStyle','--','LineWidth', 1);
xlabel('x [m]');
ylabel('Bending Moment [Nm]');
title('Bending Moment Diagram for Max-Q');
xlim([0, l_tot]);
grid on;
legend('', '', '', 'Bending Moment','Sizing Bending Moment','Location', 'eastoutside');

if SaveFlag==1
exportgraphics(gcf, 'Load Diagram.pdf', 'ContentType','vector');
end

    case 2 
STRUCT=2;
nx = FORCES.nx;
% 
%D = -nx*m*g0 + T*cos(delta);
T = (D+nx*m*g0)/cos(delta);

nz=FORCES.nz;

L = -T*sin(delta) + m*g0 + nz*g0*m;
L_cone = (0.45)*L;
L_fin = (0.55)*L;
arm = (nz*g0*m*(l_tot-x_com))/L_cone;
%arm = ((nz+1)*g0*m*(l_tot-x_com))/L_cone;

mra = l_tot-arm;
%nz = FORCES.nz;

%% Load:

%% 1.
D1 = D;

P1 = D1 + m1*g0*(nx);
R1 = -L_cone + m1*g0*nz +m1*g0;

%% 2.
P2 = P1 + m2*g0*(nx);
R2 = R1 + m2*g0*nz +m2*g0;
%% 3. 

P3 = P2 + m3*g0*(nx);
R3 = R2 + m3*g0*nz +m3*g0;
%% 4.

P4 = P3 + m4*g0*(nx);
R4 = R3 + m4*g0*nz +m4*g0;

%% 5.

P5 = P4 + m5*g0*(nx);
R5 = R4 + m5*g0*nz +m5*g0;

%% 6.

P6 = P5 + m6*g0*(nx);
R6 = R5 + m6*g0*nz +m6*g0;
%% 7.

P7 = P6 + m7*g0*(nx);
R7=R6+ m7*g0*nz +m7*g0;
%% 8.

P8 = P7 + m8*g0*(nx);
R8=R7 + m8*g0*nz +m8*g0;
%% 9.

P9 = P8 + m9*g0*(nx);
R9= R8+ m9*g0*nz +m9*g0;
%% 10.

P10 = P9 + m10*g0*(nx);
R10= R9 + m10*g0*nz+m10*g0;

P11 = P10 - T*cos(delta);
R11 = R10 - T*sin(delta)-L_fin;

% Moments:

%M_1 = 

M_cg = -L_cone*(x_com-mra) + (L_fin+T*sin(delta))*(l_tot-x_com);

%% PLOTS:

% Segment lengths
% x_b_com_v = [b1, b2, b3,b34 ,b4,b5,b6,b7];
% x_b_com = x_com-sum(x_b_com_v);
x_segments = [b1*2/3, b1/3+b2, b3 ,b4+b34,b5,b6,b7,b8+b78,b9,b10];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting

% Axial Forces (P) at segment ends
P_forces = [P1,P2,P3,P4,P5,P6,P7,P8,P10,P10];
THICKNESS.P_forces = [D,P1,P2,P3,P3,P4,P5,P6,P7,P7,P8,P9,P10];

figure;
tiledlayout(3,1)
nexttile
hold on;
for i = 1:length(x_segments)
    plot([x(i+1), x(i+1)], [D, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
end
stairs(x, [D, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Axial Force [N]');
title('Axial Force Diagram for Max-Q');
xlim([0, sum(x_segments)]);
grid on;


% Axial Forces (P) at segment ends
P_forces = -[R1,R2,R3,R4,R5,R6,R7,R8,R10,R10];
nexttile
hold on;
for i = 1:length(x_segments)
     plot([x(i+1), x(i+1)], [0, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
 end
stairs(x, [0, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Shear Force [N]');
title('Shear Force Diagram for Max-Q');
xlim([0, sum(x_segments)]);
grid on;

P_forces = [R1,R2,R3,R4,R5,R6,R7,R8,R9,R11];
x_segments = [b1*2/3, b2+b1/3, b3 ,b4+b34,b5,b6,b7,b8+b78,b9,b10];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting
    % Calculate the bending moment
M_moments = zeros(size(P_forces));
M_moments(1) = 0; % Initial bending moment

for i = 2:length(P_forces)
    M_moments(i) = M_moments(i-1) + P_forces(i-1) * x_segments(i-1);
end

M_moments(end)=0;
x_M_2 = [0,0,sum(b1+b2+b3+b4+b5+b34),sum(b1+b2+b3+b4+b5+b34),0];
y_M_2 = [0,M_moments(5),M_moments(5),0,0];
x_M_1 = [sum(b1+b2+b3+b4+b5+b34),sum(b1+b2+b3+b4+b5+b34),l_tot,l_tot,sum(b1+b2+b3+b4+b5)];
y_M_1 = [0,min(M_moments),min(M_moments),0,0];
% Plot the bending moment diagram
nexttile
hold on;
for i = 1:length(x_segments)
    plot([x(i+1), x(i+1)], [0, M_moments(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
end
plot(x,[0, M_moments], 'LineWidth', 2, 'Color', 'b'); % Step-style plot
hold on;
plot(x_M_2,y_M_2,'Color','r','LineStyle','--','LineWidth', 1);
hold on;
plot(x_M_1,y_M_1,'Color','r','LineStyle','--','LineWidth', 1);
xlabel('x [m]');
ylabel('Bending Moment [Nm]');
title('Bending Moment Diagram for Max-Q');
xlim([0, sum(x_segments)]);
grid on;
legend('', '', '','','', '', '','', '','','Bending Moment','Sizing Bending Moment','Location', 'eastoutside');
if SaveFlag==1
exportgraphics(gcf, 'Load Diagram.pdf', 'ContentType','vector');
end

    case 3

STRUCT=3;
        T = FORCES.T;

nx =((T*cos(delta) - D)/(m*g0));
%nx =((T*cos(delta) - D - m*g0)/(m*g0));

% nx = FORCES.nx;
% 
% D = -nx*m*g0 +T*cos(delta);

nz = (T*sin(delta) + L_cone + L_fin)/(m*g0);

arm = (nz*g0*m*(l_tot-x_com))/L_cone;
%arm = ((nz+1)*g0*m*(l_tot-x_com))/L_cone;

mra = l_tot-arm;
%nz = FORCES.nz;

%% Load:

%% 1.
D1 = D;

P1 = D1 + m1*g0*(nx);
R1 = -L_cone + m1*g0*nz;

%% 2.
P2 = P1 + m2*g0*(nx);
R2 = R1 + m2*g0*nz;
%% 3. 

P3 = P2 + m3*g0*(nx);
R3 = R2 + m3*g0*nz;
%% 4.

P4 = P3 + m4*g0*(nx);
R4 = R3 + m4*g0*nz-L_fin;

P5 = P4 - T*cos(delta);
R5 = R4 - T*sin(delta);

% Moments:

%M_1 = 

M_cg = -L_cone*(x_com-mra) + (L_fin+T*sin(delta))*(l_tot-x_com);

%% PLOTS:

% Segment lengths
x_segments = [b1, b2, b3,b4];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting

% Axial Forces (P) at segment ends
P_forces = [P1,P2,P3,P4];
THICKNESS.P_forces = [D,P1,P2,P3,P4];

figure;
hold on;
for i = 1:length(x_segments)
    plot([x(i+1), x(i+1)], [D, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
end
stairs(x, [D, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Axial Force [N]');
title('Axial Force Diagram for Pitch-Up Maneuver');
xlim([0, sum(x_segments)]);
grid on;

% Segment lengths
x_segments = [b1, b2, b3,b4];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting

% Axial Forces (P) at segment ends
P_forces = -[R1,R2,R3,R4];
THICKNESS.R_forces=[0,R1,R2,R3,R4];
figure;
hold on;
for i = 1:length(x_segments)
     plot([x(i+1), x(i+1)], [0, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
 end
stairs(x, [0, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Shear Force [N]');
title('Shear Force Diagram for Pitch-Up Maneuver');
xlim([0, sum(x_segments)]);
grid on;

P_forces = [R1,R2,R3,R4];
x_segments = [b1, b2, b3,b4];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting
    % Calculate the bending moment
M_moments = zeros(size(P_forces));
M_moments(1) = 0; % Initial bending moment

for i = 2:length(P_forces)
    M_moments(i) = M_moments(i-1) + P_forces(i-1) * x_segments(i-1);
end

M_moments(end)=0;

THICKNESS.moments = M_moments;

% Plot the bending moment diagram
figure;
hold on;
for i = 1:length(x_segments)
    plot([x(i+1), x(i+1)], [0, M_moments(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
end
plot(x,[0, M_moments], 'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Bending Moment [Nm]');
title('Bending Moment Diagram for Pitch-Up Maneuver');
xlim([0, sum(x_segments)]);
grid on;


end % switch
end % function