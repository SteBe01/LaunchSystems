function [THICKNESS]= Final_syst_case2nxzalpha(GEOMETRY,OUT,FORCES,MAT,AER)

% 19 pieces (no x_com)

FLAG = OUT.FLAG; % 1 max q, 2 is pitch up

b1 = GEOMETRY.b1; % fairing length
b2 = GEOMETRY.b2; % forward skirt length
b3 = GEOMETRY.b3; % height of fuel tank2
b34 = GEOMETRY.b34; % 16 cm between tanks
b4 = GEOMETRY.b4; %  height of oxygen tank2
b5 = GEOMETRY.b5; %  aft skirt length 2
b6 = GEOMETRY.b6; % interstage length
b7 = GEOMETRY.b7; % height of fuel tank1
b78 = GEOMETRY.b78; % 16 cm between tanks
b8 = GEOMETRY.b8-3; % height of oxygen tank1
b9 = GEOMETRY.b9+3; % aft skirt length 1
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
% m1 = GEOMETRY.m1/100; 
% m2 = GEOMETRY.m2/100;
% m3 = GEOMETRY.m3/100;
% m4 = GEOMETRY.m4/100;
% m5 = GEOMETRY.m5/100; 
% m6 = GEOMETRY.m6/100;
% m7 = GEOMETRY.m7/100;
% m8 = GEOMETRY.m8/100;
% m9 = GEOMETRY.m9/100;
% m10 = GEOMETRY.m10/100; 

x_com = GEOMETRY.x_com;
% x_com1 = GEOMETRY.x_com1;
% x_com2 = GEOMETRY.x_com2;
% x_com3 = GEOMETRY.x_com3;
% x_com4 = GEOMETRY.x_com4;
% x_com5 = GEOMETRY.x_com5;
% x_com6 = GEOMETRY.x_com6;
% x_com7 = GEOMETRY.x_com7;
% x_com8 = GEOMETRY.x_com8;
% x_com9 = GEOMETRY.x_com9;
% x_com10 =GEOMETRY.x_com10;

% x_fin = GEOMETRY.x_fin;
m=m1+m2+m3+m4+m5+m6+m7+m8+m9+m10;

A_ref =pi*(GEOMETRY.Diam_2^2)/4;

g0 = 9.80665; %m/s^2
FS = FORCES.FS;
q = AER.q;
Cl_nose = OUT.Cl_nose;
Cl_fin = OUT.Cl_fin;
Cd = OUT.Cd;
L_cone = q*Cl_nose*A_ref;
L_fin = q*Cl_fin*A_ref;
%D = q*Cd*A_ref;
L_cone = AER.F_nose;
L_fin = AER.F_tail/50;
D = AER.D;
W = m*g0;
alpha=OUT.alpha;

if FLAG==1

%T = FORCES.T;

nx = FORCES.nx;
% 
%D = -nx*m*g0 +T;
T = D +nx*m*g0;

%% Load:

%% 1.
D1 = D;

P1 = D1 + m1*g0*(nx);

%% 2.
P2 = P1 + m2*g0*(nx);

%% 3. 

P3 = P2 + m3*g0*(nx);

%% 4.

P4 = P3 + m4*g0*(nx);

%% 5.

P5 = P4 + m5*g0*(nx);

%% 6.

P6 = P5 + m6*g0*(nx);

%% 7.

P7 = P6 + m7*g0*(nx);

%% 8.

P8 = P7 + m8*g0*(nx);

%% 9.

P9 = P8 + m9*g0*(nx);

%% 10.

P10 = P9 + m10*g0*(nx);

P11 = P10 - T;


%%


% Segment lengths
x_segments = [b1, b2, b3 ,b4,b5,b6,b7,b8,b9,b10];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting

% Axial Forces (P) at segment ends
P_forces = [P1,P2,P3,P4,P5,P6,P7,P8,P9,P10];

figure;
hold on;
for i = 1:length(x_segments)
    plot([x(i+1), x(i+1)], [D, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
end
stairs(x, [D, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Axial Force [N]');
title('Axial Force Diagram for MAX-Q');
xlim([0, sum(x_segments)]);
grid on;

THICKNESS.FLAG=1;

elseif FLAG==2

delta = OUT.delta;

T = FORCES.T;

nx = FORCES.nx;
% 
D = -nx*m*g0 + T*cos(delta);

T = (D +nx*m*g0)/cos(delta);

nz=FORCES.nz;

L = -T*sin(delta) + m*g0 + nz*g0*m;
L_cone = (2/3)*L;
L_fin = (1/3)*L;
arm = (nz*g0*m*(l_tot-x_com))/L_cone;
%arm = ((nz+1)*g0*m*(l_tot-x_com))/L_cone;

mra = l_tot-arm;
%nz = FORCES.nz;

%% Load:

%% 1.
D1 = D;

P1 = D1 + m1*g0*(nx)+m1*g0*sin(alpha);
R1 = -L_cone + m1*g0*nz +m1*g0;

%% 2.
P2 = P1 + m2*g0*(nx)+m2*g0*sin(alpha);
R2 = R1 + m2*g0*nz +m2*g0;
%% 3. 

P3 = P2 + m3*g0*(nx)+m3*g0*sin(alpha);
R3 = R2 + m3*g0*nz +m3*g0;
%% 4.

P4 = P3 + m4*g0*(nx)+m4*g0*sin(alpha);
R4 = R3 + m4*g0*nz +m4*g0;

%% 5.

P5 = P4 + m5*g0*(nx)+m5*g0*sin(alpha);
R5 = R4 + m5*g0*nz +m5*g0;

%% 6.

P6 = P5 + m6*g0*(nx)+m6*g0*sin(alpha);
R6 = R5 + m6*g0*nz +m6*g0;
%% 7.

P7 = P6 + m7*g0*(nx)+m7*g0*sin(alpha);
R7=R6+ m7*g0*nz +m7*g0;
%% 8.

P8 = P7 + m8*g0*(nx)+m8*g0*sin(alpha);
R8=R7 + m8*g0*nz +m8*g0;
%% 9.

P9 = P8 + m9*g0*(nx)+m9*g0*sin(alpha);
R9= R8+ m9*g0*nz +m9*g0;
%% 10.

P10 = P9 + m10*g0*(nx)+m10*g0*sin(alpha);
R10= R9 + m10*g0*nz - L_fin +m10*g0;

P11 = P10 - T*cos(delta);
R11 = R10 - T*sin(delta);

% Moments:

%M_1 = 

M_cg = -L_cone*(x_com-mra) + (L_fin+T*sin(delta))*(l_tot-x_com);

%% PLOTS:

% Segment lengths
% x_b_com_v = [b1, b2, b3,b34 ,b4,b5,b6,b7];
% x_b_com = x_com-sum(x_b_com_v);
x_segments = [b1, b2, b3 ,b4,b5,b6,b7,b8,b9,b10];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting

% Axial Forces (P) at segment ends
P_forces = [P1,P2,P3,P4,P5,P6,P7,P8,P10,P10];
THICKNESS.P_forces = [D,P1,P2,P3,P3,P4,P5,P6,P7,P7,P8,P9,P10];

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


% Axial Forces (P) at segment ends
P_forces = -[R1,R2,R3,R4,R5,R6,R7,R8,R10,R10];
THICKNESS.R_forces=[0,R1,R2,R3,R3,R4,R5,R6,R7,R7,R8,R9,R10];
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

P_forces = [R1,R2,R3,R4,R5,R6,R7,R8,R9,R11];
x_segments = [b1, b2, b3 ,b4,b5,b6,b7,b8,b9,b10];
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


    THICKNESS.FLAG=2;

end % if
end