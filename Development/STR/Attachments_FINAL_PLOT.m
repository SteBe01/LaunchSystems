function [CLAMP] = Attachments_FINAL_PLOT(GEOMETRY,LOAD,PlotFlag,SaveFlag,h_it,M_it)


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
param=0.45;
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
e =d-GEOMETRY.b10;
CHECK=l_cone+z+a+b+c+d+e-L_tot;
if CHECK>10^-3

fprintf('ERROR IN LENGTH');

end
DIFF = 20;
P1 =  D + F1_x+F2_x+F3_x;
CLAMP.P1=P1;

sol = subs(sol);


if PlotFlag==1

x_segments = [l_cone+z,a,b,c,d,e];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting

% Axial Forces (P) at segment ends
P_forces = [F1_x+D,F1_x+D-P1,F1_x+D-P1+F2_x,F1_x+D-P1+F2_x,0,0];

figure()
hold on;
for i = 1:length(x_segments)
     plot([x(i+1), x(i+1)], [D, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
end
stairs(x, [D, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Axial Force [N]');
title('Axial Force Diagram for Worst Case Carrier Maneuver');
xlim([0,DIFF]);
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

figure()
hold on;
for i = 1:length(x_segments)
    plot([x(i+1), x(i+1)], [0, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
end
stairs(x, [0, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Shear Force [N]');
title('Shear Force Diagram for Worst Case Carrier Maneuver');
xlim([0,DIFF]);
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

M_moments(2)=(L_cone-(F1*2/3))*z*10^6;
M_moments(3)=0;

M_moments(4) = min(M_moments)*0.7;

M_moments(5)=0;

M_moments(end-1)=-M_moments(end-1)*15;

M_moments(end)=0;



% Plot the bending moment diagram
figure()
hold on;
for i = 1:length(x_segments)
    plot([x(i+1), x(i+1)], [0, M_moments(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
end
plot(x,[0, M_moments], 'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Bending Moment [Nm]');
title('Bending Moment Diagram for Worst Case Carrier Maneuver');
xlim([0, DIFF]);
grid on;
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
hold on;
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

end

if SaveFlag==1
exportgraphics(gcf, 'Load Diagram Attachments.pdf', 'ContentType','vector');
end

end % function