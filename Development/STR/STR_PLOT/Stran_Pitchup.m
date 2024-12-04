function [THICKNESS]= Stran_Pitchup(FORCES,GEOMETRY,MAT,OUT)

D = -(FORCES.D_x + FORCES.L_x);
L = FORCES.D_z + FORCES.L_z;

T_x = FORCES.T_x;
T_z = FORCES.T_z;

F_x = FORCES.G_x + FORCES.F_x;
F_z =FORCES.G_z + FORCES.F_z;

FS = OUT.FS;

sigma_yield_2 = MAT.MAT_2.sy;
sigma_yield_1 = MAT.MAT_1.sy;

b1 = GEOMETRY.Stk1.Cone.L;

b2 = GEOMETRY.Stk1.Stage2.L;

b3 = GEOMETRY.Stk1.Interstage.L;

b4 = GEOMETRY.Stk1.Stage1.L;

m1 =  GEOMETRY.m1;

m2 = GEOMETRY.m2;

m3 = GEOMETRY.m3;

m4 = GEOMETRY.m4;

% m1 = OUT.m/40;
% m2 = 6*OUT.m/40;
% m3 = OUT.m/40;
% m4 = 32*OUT.m/40;

Diam_1 = GEOMETRY.Diam_1;
Diam_2 = GEOMETRY.Diam_2;

l_tot = GEOMETRY.TR.Length;

g0 = OUT.g;

MAT_2 = MAT.MAT_2;
MAT_1 = MAT.MAT_1;

%% Load division:

q_D = D/l_tot;
q_L = L/l_tot;

nx = abs(F_x/(m1+m2+m3+m4));
nz =abs(F_z/(m1+m2+m3+m4));

%% 1.
D1 = D;
L1 = q_L*b1;
W1x=m1*g0*(nx);
W1z=m1*g0*(nz);

P_1 = D1 + W1x;

R_1 = L1 - W1z;

%M_1 = (L1*(b1 - x_cp_nose)) - (m1*(nz*g0 + g0*cos(alpha))*b1/3);

%% 2.
D2 = 0;
L2 = q_L*b2;
W2x=m2*g0*(nx);
W2z=m2*g0*(nz);
P_2 = P_1 + D2 +W2x;

R_2 = L2 - W2z + R_1;

%M_2 = M_1 + R_1*b2 - (m2*(g0*nz + g0*cos(alpha))*b2/2) + (L2*b2/2);

%% 3.
Rad1 = Diam_1/2; 
Rad2 = Diam_2/2;   
h = b3;    
h_cm3 = (h / 4) * ((Rad1^2 + 2*Rad1*Rad2 + 3*Rad2^2) / (Rad1^2 + Rad1*Rad2 + Rad2^2));

D3 = 0;
L3 = q_L*b3;
W3x=m3*g0*(nx);
W3z=m3*g0*(nz);
P_3 = P_2 + D3 + W3x;

R_3 = L3 -W3z + R_2;

%M_3 = M_2 + R_2*b3 - (m3*(g0*nz + g0*cos(alpha))*(b3 -h_cm3)) + (L3*b3/2);

%% 4.
D4 = 0;
L4 = q_L*b4;
W4x = m4*g0*(nx);
W4z = m4*g0*(nz);
P_4 = P_3 + D4 + W4x - T_x;

R_4 = L4 - W4z + R_3 + T_z;

%M_4 = M_3 + R_3*b4 - (m4*(g0*nz + g0*cos(alpha))*b4/2); % - (L4*b4/2)
 
%% Error computation:

M_1 = 5*10^6;
M_2 = 5*10^6;
M_3 = 5*10^6;
M_4 = 5*10^6;

Delta_M_4_err = abs(M_4);

Delta_R_4_err = (abs(R_4));

Delta_P_4_err = abs(abs(P_4));

if Delta_M_4_err~=0

fprintf('Momentum not balanced: \n Delta_M_4_err = %d [Nm] \n',abs(Delta_M_4_err));

end

if Delta_R_4_err~=0

fprintf('Shear force not balanced: \n  Delta_R_4_err = %d [N] \n',abs(Delta_R_4_err));

end

if Delta_P_4_err~=0

fprintf('Axial force not balanced: \n  Delta_P_4_err = %d [N] \n',abs(Delta_P_4_err));

end

LOAD.P_1 = P_1;
LOAD.P_2 = P_2;
LOAD.P_3 = P_3;
LOAD.P_4 = P_4;

LOAD.R_1 = R_1;
LOAD.R_2 = R_2;
LOAD.R_3 = R_3;
LOAD.R_4 = R_4;

LOAD.M_1 = M_1;
LOAD.M_2 = M_2;
LOAD.M_3 = M_3;
LOAD.M_4 = M_4;

LOAD.Delta_M_4_err = Delta_M_4_err;
LOAD.Delta_P_4_err = Delta_P_4_err;
LOAD.Delta_R_4_err = Delta_R_4_err;

%% Stresses:

% 2.

sigma_a2 = FS*abs(P_2/(2*pi*(Diam_2/2)));

sigma_b2 = FS*abs(M_2/(pi*((Diam_2/2)^2)));

sigma_max2 = sigma_a2 + sigma_b2;

%tau2 = FS*abs(R_2/A_2);

LOAD.thick_2 = sigma_max2/sigma_yield_2;

% 4.

sigma_a4 = FS*abs(P_4/(2*pi*(Diam_1/2)));

sigma_b4 = FS*abs(M_4/(pi*((Diam_1/2)^2)));

sigma_max4 = sigma_a4 + sigma_b4;

%tau4 = FS*abs(R_4/A_4);

LOAD.thick_4 = sigma_max4/sigma_yield_1;

%% Buckling:

% 2.
sigma_crit_2 = ((9*((LOAD.thick_2/(Diam_2/2))^1.6)) + (0.16*(((LOAD.thick_2/(b2))^1.3))))*MAT_2.E;

if sigma_crit_2<sigma_max2

    fprintf('Buckling failure of unpressurized 2nd stage \n');

    %t_2_n = @(x) sigma_max2 - (((9*((x/(Diam_2/2))^1.6)) + (0.16*(((x/(b2))^1.3))) + min([0.191*(p/MAT_2.E)*(((Diam_2/2)/x)^2),0.229]))*MAT_2.E*x/(Diam_2/2)); % pressurized
     t_2_n = @(x) sigma_max2 - (((9*((x/(Diam_2/2))^1.6)) + (0.16*(((x/(b2))^1.3))))*MAT_2.E*x/(Diam_2/2)); % unpressurized
    t20=LOAD.thick_2;

    LOAD.thick_2 = fsolve(t_2_n,t20);

end


% 4.
sigma_crit_1 = ((9*((LOAD.thick_4/(Diam_1/2))^1.6)) + (0.16*(((LOAD.thick_4/(b4))^1.3))))*MAT_1.E;

if sigma_crit_1<sigma_max4

    fprintf('Buckling failure of unpressurized 1st stage \n');
%t_4_n = @(x) sigma_max4 - (((9*((x/(Diam_1/2))^1.6)) + (0.16*(((x/(b4))^1.3))) + min([0.191*(p/MAT_4.E)*(((Diam_1/2)/x)^2),0.229]))*MAT_4.E*x/(Diam_1/2)); % pressurized
     t_4_n = @(x) sigma_max2 - (((9*((x/(Diam_2/2))^1.6)) + (0.16*(((x/(b2))^1.3))))*MAT_2.E*x/(Diam_2/2)); % unpressurized
 t40=LOAD.thick_4;
LOAD.thick_4 = fsolve(t_4_n,t40);
end

% 3.

alpha_c = atan(((Diam_2/2)-(Diam_1/2))/(b3));

P_cr = 0.33*(2*pi*MAT_2.E*(LOAD.thick_2^2)*(cos(alpha_c)^2))/(sqrt(3*(1-MAT_2.nu^2)));

M_cr = 0.41*(2*pi*MAT_2.E*(LOAD.thick_2^2)*(cos(alpha_c)^2)*(Diam_2/2))/(sqrt(3*(1-MAT_2.nu^2)));

LOAD.thick_3 = LOAD.thick_2;

if P_cr<P_3 && M_cr<M_3

fprintf('Buckling failure of unpressurized Interstage \n');    

P_cr = 0.33*(2*pi*MAT_2.E*(LOAD.thick_4^2)*(cos(alpha_c)^2))/(sqrt(3*(1-MAT_2.nu^2)));

M_cr = 0.41*(2*pi*MAT_2.E*(LOAD.thick_4^2)*(cos(alpha_c)^2)*(Diam_2/2))/(sqrt(3*(1-MAT_2.nu^2)));

if P_cr<P_3 && M_cr<M_3

t_3_n = @(x) ((P_3/(0.33*(2*pi*MAT_2.E*(x^2)*(cos(alpha_c)^2))/(sqrt(3*(1-MAT_2.nu^2))))) + (M_3/(0.41*(2*pi*MAT_2.E*(x^2)*(cos(alpha_c)^2)*(Diam_2/2))/(sqrt(3*(1-MAT_2.nu^2)))))) - 1;

t30=LOAD.thick_4;
LOAD.thick_3 = fsolve(t_3_n,t30);

end

end

LOAD.thick_1 =  LOAD.thick_2;

if LOAD.thick_2<MAT_2.t_min

LOAD.thick_2 = MAT_2.t_min;
LOAD.thick_1 = MAT_2.t_min;


end

if LOAD.thick_3<MAT_1.t_min

LOAD.thick_3 = MAT_1.t_min; 


end

if LOAD.thick_4<MAT_1.t_min

LOAD.thick_4 = MAT_1.t_min; 


end

D_int_base = Diam_2 - 2*LOAD.thick_2;
V_int_frustum = pi*(b3/3)*(Rad2^2 + Rad2*Rad1 + Rad1^2);
V_ext_frustum = pi*(b3/3)*((Rad2+LOAD.thick_3)^2 + (Rad2+LOAD.thick_3)*(Rad1+LOAD.thick_3) + (Rad1+LOAD.thick_3)^2);
V_fair = V_fair_approx_thickness(D_int_base,GEOMETRY.TR.Fair.fair_length*2/5,D_int_base/2,GEOMETRY.TR.Fair.fair_length*2/5,GEOMETRY.TR.Fair.fair_length/5,LOAD.thick_2);
LOAD.m1_new = (MAT_2.rho)*V_fair;
LOAD.m2_new = (MAT_2.rho)*b2*pi*Diam_2*LOAD.thick_2;
LOAD.m3_new = (MAT_2.rho)*(V_ext_frustum - V_int_frustum);
LOAD.m4_new = (MAT_1.rho)*b4*pi*Diam_1*LOAD.thick_4;

LOAD.Ms2 = LOAD.m1_new + LOAD.m2_new;
LOAD.Ms1 = LOAD.m3_new + LOAD.m4_new;

%% PLOT:
if nargin>1
% Step Plot for Axial Forces (P)

% Segment lengths
x_segments = [b1, b2, b3, b4];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting

% Axial Forces (P) at segment ends
P_forces = [P_2, P_3, P_4, P_4];

figure;
hold on;
for i = 1:length(x_segments)
    plot([x(i+1), x(i+1)], [0, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
end
stairs(x, [P_1, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Axial Force [N]');
title('Axial Force Diagram for Pitch-up Maneuver');
xlim([0, b1+b2+b3+b4]);
grid on;

% Step Plot for Shear Forces (R)
% Segment lengths
x_segments = [b1, b2, b3, b4];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plottings
P_forces = [R_2, R_3, R_4, R_4];
figure;
hold on;
for i = 1:length(x_segments)
    plot([x(i+1), x(i+1)], [0, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
end
stairs(x, [R_1, P_forces], 'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Shear Force [N]');
title('Shear Force Diagram for Pitch-up Maneuver');
xlim([0, b1+b2+b3+b4]);
grid on;

% Segment lengths
x_segments = [b1, b2, b3, b4]; % Lengths of individual segments
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting

% Moments (M) at segment ends
M_forces = [M_1, M_2, M_3, M_4]; % Moments at the end of each segment

% Plot
figure;
plot(x, [0, M_forces], 'LineWidth', 2, 'Color', 'b'); % Include zero at start
xlabel('x [m]');
ylabel('Moment [Nm]');
title('Moment Diagram for Pitch-up Maneuver');
xlim([0, b1+b2+b3+b4]);
grid on;

% Add dashed vertical lines for segment boundaries
hold on;
for i = 1:length(x_segments)
    plot([x(i+1), x(i+1)], [0, M_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
end
end


end