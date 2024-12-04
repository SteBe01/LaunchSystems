function [THICKNESS]= Final_syst_case3(GEOMETRY,OUT,FORCES,MAT,AER)

% 4 pieces

FLAG = OUT.FLAG; % 1 max q, 2 is pitch up

b1 = GEOMETRY.b1; % fairing length
b2 = GEOMETRY.b2+ GEOMETRY.b3+ GEOMETRY.b3+GEOMETRY.b4+GEOMETRY.b5; % forward skirt length
b3 = GEOMETRY.b6; % interstage length
b4 = GEOMETRY.b7+GEOMETRY.b78+GEOMETRY.b8+GEOMETRY.b9+GEOMETRY.b10; % height of fuel tank1

x_segments = [b1, b2, b3,b4];
l_tot=sum(x_segments);

m1 = GEOMETRY.m1; 
m2 = GEOMETRY.m2+GEOMETRY.m3+GEOMETRY.m4+GEOMETRY.m5;
m3 = GEOMETRY.m6;
m4 = GEOMETRY.m7+GEOMETRY.m8+GEOMETRY.m9+GEOMETRY.m10;

x_com = GEOMETRY.x_com;

% x_fin = GEOMETRY.x_fin;
m=m1+m2+m3+m4;

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
L_fin = AER.F_tail;
D = AER.D;
W = m*g0;

if FLAG==1

T = FORCES.T;

nx =((T - D)/(m*g0));

%nx =((T - D-m*g0)/(m*g0));

% nx = FORCES.nx;
% 
% D = -nx*m*g0 +T;

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


P5 = P4- T;


%%


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
title('Axial Force Diagram for MAX-Q');
xlim([0, sum(x_segments)]);
grid on;

THICKNESS.FLAG=1;

elseif FLAG==2

delta = OUT.delta;

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


    THICKNESS.FLAG=2;

end % if
end