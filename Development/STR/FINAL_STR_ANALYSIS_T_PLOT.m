function [STRUCT]=FINAL_STR_ANALYSIS_T_PLOT(GEOMETRY,FORCES,CASE,SaveFlag,h_it,M_it)

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

%[x_com] = MASS_PROPERTIES_1(GEOMETRY.m_prop_1,GEOMETRY);
x_com = GEOMETRY.x_com;
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
x_segments = [x_com,l_tot-x_com-b10];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting

DIFF = 20;

% Axial Forces (P) at segment ends
P_forces = [P2,P2];

figure();
for i = 1:length(x_segments)
    plot([x(i+1), x(i+1)], [D, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
end
stairs(x, [D, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Axial Force [N]');
title('Axial Force Diagram for Max-Q');
xlim([0, DIFF]);
grid on;

% Segment lengths
x_segments = [b1*2/3,x_com-(b1*2/3),l_tot-x_com-b10];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting

% Axial Forces (P) at segment ends
P_forces = -[R1,R2,R2];
figure()
for i = 1:length(x_segments)
     plot([x(i+1), x(i+1)], [0, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
 end
stairs(x, [0, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Shear Force [N]');
title('Shear Force Diagram for Max-Q');
xlim([0, DIFF]);
grid on;

P_forces = [R1,R2,R2];
x_segments =[b1*2/3,x_com-(b1*2/3),l_tot-x_com-b10] ;
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
x_M_1 = [sum(b1+b2+b3+b4+b5+b34),sum(b1+b2+b3+b4+b5+b34),l_tot-b10,l_tot-b10,sum(b1+b2+b3+b4+b5)];
y_M_1 = [0,min(M_moments),min(M_moments),0,0];

% Plot the bending moment diagram
figure()
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
xlim([0, DIFF]);
grid on;
legend('', '', '', 'Bending Moment','Sizing Bending Moment','Location', 'eastoutside');
%yticks = get(gca, 'YTick');                         
% yticklabels = arrayfun(@(v) sprintf('%d \\times 10^5', abs(v) / 1e5), yticks, 'UniformOutput', false);
% set(gca, 'YTickLabel', yticklabels);
yticks = get(gca, 'YTick'); % Get current tick values
% Convert tick values to scientific notation strings with "e" notation and absolute values
yticklabels = arrayfun(@(y) sprintf('%0.2e', abs(y)), yticks, 'UniformOutput', false);
% Set the formatted tick labels to the y-axis
set(gca, 'YTickLabel', yticklabels);
%fundamental parameters recovery
AR = h_it.stg1.R_lox / h_it.stg1.dome_lox; %[-] aspect ratio of the tanks
diam1 = M_it.diam1; %[m] first stage diameter
r1 = diam1/2; %[m] first stage radius
diam2 = M_it.diam2; %[m] second stage diameter
r2 = diam2/2; %[m] second stage radius

%height recovery
%first stage:
h0 = 0; %[m] starting height 
h1 = h_it.stg1.motor; %[m] height of the first stage motor
h2 = h_it.stg1.C3 + h1; %[m] height at which the cylindrical part of the lox tank starts (first stage)
h3 = h2 + h_it.stg1.cyl_lox +...
    h_it.stg1.C2 + h_it.stg1.cyl_rp1; %[m] height at which the cylindrical part of the rp1 tank ends (first stage)
%h4 = h_it(r).stg1.tot; %[m] height of the first stage end
h4 = h3 + h_it.stg1.C1; %[m] height of the first stage end
%second stage:
h5 = h4 + h_it.stg2.C3; %[m] height at which the cylindrical part of the lox tank starts (second stage)
h6 = h5 + h_it.stg2.cyl_lox +...
    h_it.stg2.C2 + h_it.stg2.cyl_rp1; %[m] height at which the cylindrical part of the rp1 tank ends (second stage)
h7 = h6 + h_it.stg2.C1; %[m] height of the second stage end
%fairing:
% h8 = h7 + h_it(r).fairing - 2 * diam2; %[m] height at which the fairing cone starts
h8 = h7; %[m] height at which the fairing cone starts
h9 = h7 + h_it.fairing; %[m] total height of the rocket
%CoM:
hCG = h_it.CG; %[m] total rocket CoM position (in "h" coordinates)

%silhouette parameters:
h_rocket = [h1, h1, h2, h3, h4, h5, h6, h7, h8, h9]';
y = [ 0, r1, r1, r1, r2, r2, r2, r2, r2,  0]'; %[m]
h_rocket = [h_rocket; flip(h_rocket)]; %[m] change coordinate system to adequate to the convention
x = h9 - h_rocket; %[m]
y = [y; -flip(y)]; %[m]
xCG = h9 - hCG; %[m] total rocket CoM position (in "x" coordinates)
yCG = 0; %[m]

%silhouette plot:
nexttile;
hold on;
plot(x, y, '-k', LineWidth=1); grid on; axis equal;hold on;
% plot([diam1/2, -diam1/2], [h1.til_tank-h1.dome_rp1,h1.til_tank-h1.dome_rp1], '--k');
% plot([diam2/2, -diam2/2], [h1.attach,h1.attach], '--k');
% plot([diam2/2, -diam2/2], [h3.tot-2*diam2,h3.tot-2*diam2], '--k');
plot(xCG, yCG, '+r'); 

%tank plot:
cap1 = @(k) sqrt(r1^2 - k.^2)/AR;
y1 = linspace(-r1, r1, 1e4);
x1 = @(k) h9 - ( h2 - cap1(k) );
plot(x1(y1), y1, '-k');
y2 = y1;
x2 = @(k) h9 - ( h2 + h_it.stg1.cyl_lox + cap1(k) );
plot(x2(y2), y2, '-k');
AR_1 = ( r1 + h_it.stg1.C2 ) / ( h_it.stg1.dome_lox + h_it.stg1.C2);
cap1_AR_1 = @(k) sqrt( ( r1 + h_it.stg1.C2 )^2 - k.^2)/AR_1;
y3 = y1;%linspace(-r1 - h_it(r).stg1.C2, r1 + h_it(r).stg1.C2, 1e4); 
x3 = @(k) h9 - ( h2 + h_it.stg1.cyl_lox + cap1_AR_1(k) );
plot(x3(y3), y3, '-k');
y4 = y1;
x4 = @(k) h9 - ( h3 + cap1(k) );
plot(x4(y4), y4, '-k');
cap2 = @(k) sqrt(r2^2 - k.^2)/AR;
y5 = linspace(-r2, r2, 1e4);
x5 = @(k) h9 - ( h5 - cap2(k) );
plot(x5(y5), y5, '-k');
y6 = y5;
x6 = @(k) h9 - ( h5 + h_it.stg2.cyl_lox + cap2(k) );
plot(x6(y6), y6, '-k');
AR_2 = ( r2 + h_it.stg2.C2 ) / ( h_it.stg2.dome_lox + h_it.stg2.C2);
cap2_AR_2 = @(k) sqrt( ( r2 + h_it.stg2.C2 )^2 - k.^2)/AR_2;
y7 = y5;%linspace(-r2 - h_it(r).stg2.C2, r2 + h_it(r).stg2.C2, 1e4); 
x7 = @(k) h9 - ( h5 + h_it.stg2.cyl_lox + cap2_AR_2(k) );
plot(x7(y7), y7, '-k');
y8 = y5;
x8 = @(k) h9 - ( h6 + cap2(k) );
plot(x8(y8), y8, '-k');

%fairing plot:
x9 = h9 - h8 * ones(2,1);
y9 = [-r2; r2];
plot(x9, y9, '--k',LineWidth=1);
x10 = h9 - h7 * ones(2,1);
y10 = y9;
plot(x10, y10, '--k',LineWidth=1);

%motor #i:
xm = [0, 0, (h1-0.26)*0.5, h1-0.26, h1-0.19, h1, h1, h1-0.19, h1-0.26, (h1-0.26)*0.5, 0, 0];
h1_2 = h_it.stg2.motor;
xm_2 = [0, 0, (h1_2-0.26)*0.5, h1_2-0.26, h1_2-0.19, h1_2, h1_2, h1_2-0.19, h1_2-0.26, (h1_2-0.26)*0.5, 0, 0];
ym = [0, 0.15, 0.12, 0.035, 0.063, 0.063, -0.063, -0.063, -0.035, -0.12, -0.15, 0];
grayColor = [128 128 128]/255;

%first stage:
alpha = linspace(0, 2*pi, M_it.n_mot1 + 1);
for i = 1 : M_it.n_mot1
    xm1 = h9 - xm;
    ym1 = (r1 - 0.15)*cos(alpha(i)) - ym;
    plot(xm1, ym1, '-','Color', grayColor, LineWidth=1);
end

%second stage:
xm2 = h9 - ( h4 + xm_2 - h1_2 );
ym2 = ym;
plot(xm2, ym2, '-','Color', grayColor, LineWidth=1);
ystg2 = y9;
xstg2 = h9 - h_it.stg1.tot * ones(2,1);
plot(xstg2, ystg2, '-k',LineWidth=1);
 
%adapter:
h_ad = r2 / 2;
r1_ad = r2 / 1.6;
r2_ad = r2 / 1.2;
X_a = (h9 - h7) - [0, 0, h_ad, h_ad, 0, 0];
Y_a = [0, r2_ad, r1_ad, -r1_ad, -r2_ad, 0];
plot(X_a, Y_a, '--r'); %adapter

%payload
r_pl = 0.90 * r1_ad;
theta = linspace(0, 2*pi, 3e1 + 1);
ypl = r_pl * sin(theta);
xpl = (h9 - h7 - h_ad * 0.3 - r_pl ) + r_pl * cos(theta);
plot( xpl, ypl, '--b');
xlim([0;DIFF])

if SaveFlag==1
exportgraphics(gcf,'Load Diagram.pdf', 'ContentType','vector');
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

P10 = P9 + m10*g0*(nx)-T*cos(delta);
R10= R9 + m10*g0*nz+m10*g0- T*sin(delta)-L_fin;

% P11 = P10 - T*cos(delta);
% R11 = R10 - T*sin(delta)-L_fin;

% Moments:

%M_1 = 

M_cg = -L_cone*(x_com-mra) + (L_fin+T*sin(delta))*(l_tot-x_com);

%% PLOTS:

% Segment lengths
% x_b_com_v = [b1, b2, b3,b34 ,b4,b5,b6,b7];
% x_b_com = x_com-sum(x_b_com_v);
x_segments = [b1*2/3, b1/3+b2, b3 ,b4+b34,b5,b6,b7,b8+b78,b9];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting
DIFF =20;
% Axial Forces (P) at segment ends
P_forces = [P1,P2,P3,P4,P5,P6,P7,P9,P9];

figure()
hold on;
for i = 1:length(x_segments)
    plot([x(i+1), x(i+1)], [D, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
end
stairs(x, [D, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Axial Force [N]');
title('Axial Force Diagram for Max-Q');
xlim([0, DIFF]);
grid on;


% Axial Forces (P) at segment ends
P_forces = -[R1,R2,R3,R4,R5,R6,R7,R9,R9];
figure()
hold on;
for i = 1:length(x_segments)
     plot([x(i+1), x(i+1)], [0, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
 end
stairs(x, [0, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Shear Force [N]');
title('Shear Force Diagram for Max-Q');
xlim([0, DIFF]);
grid on;

P_forces = [R1,R2,R3,R4,R5,R6,R7,R8,R9];
x_segments = [b1*2/3, b2+b1/3, b3 ,b4+b34,b5,b6,b7,b8+b78,b9];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting
    % Calculate the bending moment
M_moments = zeros(size(P_forces));
M_moments(1) = 0; % Initial bending moment

for i = 2:length(P_forces)
    M_moments(i) = M_moments(i-1) + P_forces(i-1) * x_segments(i-1);
end

M_interp_2 = [0;M_moments(3)];
x_interp_2 = [b1*2/3;sum(b1+b2+b3)];
M_fit =polyfit(x_interp_2,M_interp_2,1);
M_moments(2)= M_fit(1)*(b1+b2)+M_fit(2);

M_moments(end-1)=M_moments(end-1)*1;
M_moments(end)=0;
x_M_2 = [0,0,sum(b1+b2+b3+b4+b5+b34),sum(b1+b2+b3+b4+b5+b34),0];
y_M_2 = [0,M_moments(5),M_moments(5),0,0];
x_M_1 = [sum(b1+b2+b3+b4+b5+b34),sum(b1+b2+b3+b4+b5+b34),l_tot-b10,l_tot-b10,sum(b1+b2+b3+b4+b5)];
y_M_1 = [0,min(M_moments),min(M_moments),0,0];
% Plot the bending moment diagram
figure()
hold on;
for i = 1:length(x_segments)
    plot([x(i+1), x(i+1)], [0, M_moments(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
end
plot(x,[0, M_moments], 'LineWidth', 2, 'Color', 'b'); % Step-style plot
hold on;
plot([0;sum(b1+b2+b3+b4+b5+b34)],[0;0],'Color','r','LineStyle','--','LineWidth', 1);
hold on;
plot([0;0],[0;M_moments(5)],'Color','r','LineStyle','--','LineWidth', 1);
hold on;
plot([0;sum(b1+b2+b3+b4+b5+b34)],[M_moments(5);M_moments(5)],'Color','r','LineStyle','--','LineWidth', 1);
hold on;
plot(x_M_1,y_M_1,'Color','r','LineStyle','--','LineWidth', 1);
xlabel('x [m]');
ylabel('Bending Moment [Nm]');
title('Bending Moment Diagram for Max-Q');
xlim([0, DIFF]);
grid on;
legend('', '', '','','', '', '','', '','Bending Moment','Sizing Bending Moment','Location', 'best');
%yticks = get(gca, 'YTick');                         
% yticklabels = arrayfun(@(v) sprintf('%d \\times 10^5', abs(v) / 1e5), yticks, 'UniformOutput', false);
% set(gca, 'YTickLabel', yticklabels);
yticks = get(gca, 'YTick'); % Get current tick values

% Convert tick values to scientific notation strings with "e" notation and absolute values
yticklabels = arrayfun(@(y) sprintf('%0.2e', abs(y)), yticks, 'UniformOutput', false);

% Set the formatted tick labels to the y-axis
set(gca, 'YTickLabel', yticklabels);
%fundamental parameters recovery
AR = h_it.stg1.R_lox / h_it.stg1.dome_lox; %[-] aspect ratio of the tanks
diam1 = M_it.diam1; %[m] first stage diameter
r1 = diam1/2; %[m] first stage radius
diam2 = M_it.diam2; %[m] second stage diameter
r2 = diam2/2; %[m] second stage radius

%height recovery
%first stage:
h0 = 0; %[m] starting height 
h1 = h_it.stg1.motor; %[m] height of the first stage motor
h2 = h_it.stg1.C3 + h1; %[m] height at which the cylindrical part of the lox tank starts (first stage)
h3 = h2 + h_it.stg1.cyl_lox +...
    h_it.stg1.C2 + h_it.stg1.cyl_rp1; %[m] height at which the cylindrical part of the rp1 tank ends (first stage)
%h4 = h_it(r).stg1.tot; %[m] height of the first stage end
h4 = h3 + h_it.stg1.C1; %[m] height of the first stage end
%second stage:
h5 = h4 + h_it.stg2.C3; %[m] height at which the cylindrical part of the lox tank starts (second stage)
h6 = h5 + h_it.stg2.cyl_lox +...
    h_it.stg2.C2 + h_it.stg2.cyl_rp1; %[m] height at which the cylindrical part of the rp1 tank ends (second stage)
h7 = h6 + h_it.stg2.C1; %[m] height of the second stage end
%fairing:
% h8 = h7 + h_it(r).fairing - 2 * diam2; %[m] height at which the fairing cone starts
h8 = h7; %[m] height at which the fairing cone starts
h9 = h7 + h_it.fairing; %[m] total height of the rocket
%CoM:
hCG = h_it.CG; %[m] total rocket CoM position (in "h" coordinates)

%silhouette parameters:
h_rocket = [h1, h1, h2, h3, h4, h5, h6, h7, h8, h9]';
y = [ 0, r1, r1, r1, r2, r2, r2, r2, r2,  0]'; %[m]
h_rocket = [h_rocket; flip(h_rocket)]; %[m] change coordinate system to adequate to the convention
x = h9 - h_rocket; %[m]
y = [y; -flip(y)]; %[m]
xCG = h9 - hCG; %[m] total rocket CoM position (in "x" coordinates)
yCG = 0; %[m]

%silhouette plot:
figure()
plot(x, y, '-k', LineWidth=1, HandleVisibility='off'); grid on; axis equal;hold on;
% plot([diam1/2, -diam1/2], [h1.til_tank-h1.dome_rp1,h1.til_tank-h1.dome_rp1], '--k');
% plot([diam2/2, -diam2/2], [h1.attach,h1.attach], '--k');
% plot([diam2/2, -diam2/2], [h3.tot-2*diam2,h3.tot-2*diam2], '--k');
plot(xCG, yCG, '+r'); 

%tank plot:
cap1 = @(k) sqrt(r1^2 - k.^2)/AR;
y1 = linspace(-r1, r1, 1e4);
x1 = @(k) h9 - ( h2 - cap1(k) );
plot(x1(y1), y1, '-k', HandleVisibility='off');
y2 = y1;
x2 = @(k) h9 - ( h2 + h_it.stg1.cyl_lox + cap1(k) );
plot(x2(y2), y2, '-k', HandleVisibility='off');
AR_1 = ( r1 + h_it.stg1.C2 ) / ( h_it.stg1.dome_lox + h_it.stg1.C2);
cap1_AR_1 = @(k) sqrt( ( r1 + h_it.stg1.C2 )^2 - k.^2)/AR_1;
y3 = y1;%linspace(-r1 - h_it(r).stg1.C2, r1 + h_it(r).stg1.C2, 1e4); 
x3 = @(k) h9 - ( h2 + h_it.stg1.cyl_lox + cap1_AR_1(k) );
plot(x3(y3), y3, '-k', HandleVisibility='off');
y4 = y1;
x4 = @(k) h9 - ( h3 + cap1(k) );
plot(x4(y4), y4, '-k');
cap2 = @(k) sqrt(r2^2 - k.^2)/AR;
y5 = linspace(-r2, r2, 1e4);
x5 = @(k) h9 - ( h5 - cap2(k) );
plot(x5(y5), y5, '-k');
y6 = y5;
x6 = @(k) h9 - ( h5 + h_it.stg2.cyl_lox + cap2(k) );
plot(x6(y6), y6, '-k', HandleVisibility='off');
AR_2 = ( r2 + h_it.stg2.C2 ) / ( h_it.stg2.dome_lox + h_it.stg2.C2);
cap2_AR_2 = @(k) sqrt( ( r2 + h_it.stg2.C2 )^2 - k.^2)/AR_2;
y7 = y5;%linspace(-r2 - h_it(r).stg2.C2, r2 + h_it(r).stg2.C2, 1e4); 
x7 = @(k) h9 - ( h5 + h_it.stg2.cyl_lox + cap2_AR_2(k) );
plot(x7(y7), y7, '-k', HandleVisibility='off');
y8 = y5;
x8 = @(k) h9 - ( h6 + cap2(k) );
plot(x8(y8), y8, '-k', HandleVisibility='off');

%fairing plot:
x9 = h9 - h8 * ones(2,1);
y9 = [-r2; r2];
plot(x9, y9, '--k',LineWidth=1, HandleVisibility='off');
x10 = h9 - h7 * ones(2,1);
y10 = y9;
plot(x10, y10, '--k',LineWidth=1, HandleVisibility='off');

%motor #i:
xm = [0, 0, (h1-0.26)*0.5, h1-0.26, h1-0.19, h1, h1, h1-0.19, h1-0.26, (h1-0.26)*0.5, 0, 0];
h1_2 = h_it.stg2.motor;
xm_2 = [0, 0, (h1_2-0.26)*0.5, h1_2-0.26, h1_2-0.19, h1_2, h1_2, h1_2-0.19, h1_2-0.26, (h1_2-0.26)*0.5, 0, 0];
ym = [0, 0.15, 0.12, 0.035, 0.063, 0.063, -0.063, -0.063, -0.035, -0.12, -0.15, 0];
grayColor = [128 128 128]/255;

%first stage:
alpha = linspace(0, 2*pi, M_it.n_mot1 + 1);
for i = 1 : M_it.n_mot1
    xm1 = h9 - xm;
    ym1 = (r1 - 0.15)*cos(alpha(i)) - ym;
    plot(xm1, ym1, '-','Color', grayColor, LineWidth=1, HandleVisibility='off');
end

%second stage:
xm2 = h9 - ( h4 + xm_2 - h1_2 );
ym2 = ym;
plot(xm2, ym2, '-','Color', grayColor, LineWidth=1, HandleVisibility='off');
ystg2 = y9;
xstg2 = h9 - h_it.stg1.tot * ones(2,1);
plot(xstg2, ystg2, '-k',LineWidth=1, HandleVisibility='off');
 
%adapter:
h_ad = r2 / 2;
r1_ad = r2 / 1.6;
r2_ad = r2 / 1.2;
X_a = (h9 - h7) - [0, 0, h_ad, h_ad, 0, 0];
Y_a = [0, r2_ad, r1_ad, -r1_ad, -r2_ad, 0];
plot(X_a, Y_a, '--r', HandleVisibility='off'); %adapter

%payload
r_pl = 0.90 * r1_ad;
theta = linspace(0, 2*pi, 3e1 + 1);
ypl = r_pl * sin(theta);
xpl = (h9 - h7 - h_ad * 0.3 - r_pl ) + r_pl * cos(theta);
plot( xpl, ypl, '--b', HandleVisibility='off');
xlim([0;DIFF])
legend('Global C.O.M.');


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