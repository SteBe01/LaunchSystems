function [LOAD,TR] = Max_Q_analysis(COM,TR,loads,LOAD,PlotFlag)

%% Definition of distributed aerodynamical loads as uniform (rectangular) distribution

g0 = 9.80665; %m/s^2

mat_id_1 = 2;
mat_id_2 = 4;

[MAT_1] = material_selection(mat_id_1);
[MAT_2] = material_selection(mat_id_2);

pressure_type = 2;

switch pressure_type 
    case 0 % case for unpressurized vessel
        p = 0;
    case 1 %case for pressure-fed
        p = 2.8 * 1e6; %[Pa] internal tank pressure
    case 2 %case for pump-fed
        p = 300 * 1e3; %[Pa] internal tank pressure
    case 3 % case for blowdown
        p = 50 * 1e6; %[Pa] internal tank pressure
end



FS = LOAD.FS;
sigma_yield_2 = MAT_2.sy;
sigma_yield_1 = MAT_1.sy;

q = 0.5*loads.rho*loads.v^2;

q = loads.q;

A_ref = (TR.Diameter*TR.Length)*pi;

D = q*A_ref*loads.Cd;

q_D = D/TR.Length;

nx = loads.nx;

T = loads.T;

%% Definition of geometry:

b1 = LOAD.Stk1.Cone.L;

b2 = LOAD.Stk1.Stage2.L;

b3 = LOAD.Stk1.Interstage.L;

b4 = LOAD.Stk1.Stage1.L;

m1 =  LOAD.Stk1.Cone.m_cone;

m2 = LOAD.Stk1.Stage2.m_stage2;

m3 = LOAD.Stk1.Interstage.m_interstage;

m4 = LOAD.Stk1.Stage1.m_stage1;

Diam_1 = TR.Diam_1;
Diam_2 = TR.Diam_2;

l_tot = TR.Length;

x_com_tot = COM.Stk1.TR.x_LV;
x_com_tot0 = COM.Stk1.TR.x_LV0;
m_tot = COM.Stk1.TR.M_LV;
m_tot0 = COM.Stk1.TR.M_LV0;

%% Load analysis:

% Verify that global balance is correct: 

P_balance = T - D - m_tot*g0*(nx + 1);

if P_balance> 1

nx = (T-D)/(m_tot*g0) - 1;

end

%% 1.
D1 = q_D*b1;

P_1 = D1 + m1*g0*(1+nx);


%% 2.
D2 = q_D*b2;

P_2 = P_1 + D2 + m2*g0*(nx + 1);

%% 3.
Rad1 = Diam_1/2; 
Rad2 = Diam_2/2;   
D3 = q_D*b3;

P_3 = P_2 + D3 + m3*g0*(nx+ 1);

%% 4.
D4 = q_D*b4;

P_4 = P_3 + D4 + m4*g0*(nx+ 1) - T;

%% Error computation:


Delta_P_4_err = abs(abs(P_4));

if Delta_P_4_err~=0

fprintf('Axial force not balanced: \n  Delta_P_4_err = %d [N] \n',abs(Delta_P_4_err));

end

LOAD.P_1 = P_1;
LOAD.P_2 = P_2;
LOAD.P_3 = P_3;
LOAD.P_4 = P_4;

LOAD.Delta_P_4_err = Delta_P_4_err;

%% Stresses:

P_1 = 4*10^5;
P_2 = 4*10^5;
P_3 = 4*10^5;
P_4 = 4*10^5;

% 2.

sigma_a2 = FS*abs(P_2/(2*pi*(Diam_2/2)));

sigma_max2 = sigma_a2;

LOAD.thick_2 = sigma_max2/sigma_yield_2;

% 4.

sigma_a4 = FS*abs(P_4/(2*pi*(Diam_1/2)));

sigma_max4 = sigma_a4;

LOAD.thick_4 = sigma_max4/sigma_yield_1;

%% Buckling:

% 2.

P_cr_2 = ((0.5625*MAT_2.E)/(4*(1-(MAT_2.nu^2))))*((LOAD.thick_2/Rad2)^3);

if P_cr_2<P_2 

fprintf('Buckling failure of unpressurized Stage 2 \n');    

LOAD.thick_2 = ((((P_cr_2*(4*(1-(MAT_2.nu^2))))/(0.5625*MAT_2.E)))^(1/3))*Rad2;

end

% 3.

alpha_c = atan(((Diam_2/2)-(Diam_1/2))/(b3));

P_cr_3 = 0.33*(2*pi*MAT_2.E*(LOAD.thick_2^2)*(cos(alpha_c)^2))/(sqrt(3*(1-MAT_2.nu^2)));

LOAD.thick_3 = LOAD.thick_2;

if P_cr_3<P_3 

fprintf('Buckling failure of unpressurized Interstage \n');    


LOAD.thick_3 = sqrt(((sqrt(3*(1-MAT_2.nu^2)))*P_cr_3)/(0.33*2*pi*MAT_2.E*(cos(alpha_c)^2)));

end

% 4.
P_cr_4 = ((0.5625*MAT_1.E)/(4*(1-(MAT_1.nu^2))))*((LOAD.thick_4/Rad1)^3);

if P_cr_4<P_4 

fprintf('Buckling failure of unpressurized Stage 1 \n');    

LOAD.thick_4 = ((((P_cr_4*(4*(1-(MAT_1.nu^2))))/(0.5625*MAT_1.E)))^(1/3))*Rad1;

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
V_fair = V_fair_approx_thickness(D_int_base,TR.Fair.fair_length*2/5,D_int_base/2,TR.Fair.fair_length*2/5,TR.Fair.fair_length/5,LOAD.thick_2);
LOAD.m1_new = (MAT_2.rho)*V_fair;
LOAD.m2_new = (MAT_2.rho)*b2*pi*Diam_2*LOAD.thick_2;
LOAD.m3_new = (MAT_2.rho)*(V_ext_frustum - V_int_frustum);
LOAD.m4_new = (MAT_1.rho)*b4*pi*Diam_1*LOAD.thick_4;

LOAD.Ms2 = LOAD.m1_new + LOAD.m2_new;
LOAD.Ms1 = LOAD.m3_new + LOAD.m4_new;

%% PLOT:
if PlotFlag==1
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
title('Axial Force Diagram for MAX-Q');
xlim([0, b1+b2+b3+b4]);
grid on;

end % function