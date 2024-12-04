function [BALANCE]= Load_Distr_BODY(OUT,GEOMETRY)

% switch CASE
% 
%     case 2 % Pitch-up 

X_balance = OUT.D + OUT.T_x - OUT.F_x + OUT.W_x;
Z_balance = OUT.L + OUT.T_z - OUT.F_z + OUT.W_z;

D=OUT.D;
L=OUT.L;
Tx=OUT.T_x;
Tz=OUT.T_z;
Fx=OUT.F_x;
Fz=OUT.F_z;
Wx=OUT.W_x;
Wz = OUT.W_z;

M=OUT.moment;
M4 =M;
T_M=OUT.T_M;
D_M=OUT.D_M;
L_M=OUT.L_M;
Aer_arm=OUT.Aer_arm;
T_arm=OUT.T_arm;
A_M=OUT.A_M; % D*sin(alpha)-L*cos(alpha)

M_balance = M -T_M*T_arm -D_M*Aer_arm +L_M*Aer_arm;
M_balance2 = M -T_M*T_arm -A_M*Aer_arm;

 if X_balance>10^-4
printf('Error in X Balance');
 end

  if Z_balance>10^-4
printf('Error in Z Balance');
  end

   if  M_balance>10^-4 && M_balance2>10^-4
printf('Error in M Balance');
 end

b1 = GEOMETRY.b1; % x_c
b2 = GEOMETRY.b2;
b3 = GEOMETRY.b3;
b4 = GEOMETRY.b4;
m1 = GEOMETRY.m1;
m2 = GEOMETRY.m2;
m3 = GEOMETRY.m3;
m4 = GEOMETRY.m4;
x_com = GEOMETRY.x_com;
x_com1 = GEOMETRY.x_com1;
x_com2 = GEOMETRY.x_com2;
x_com3 = GEOMETRY.x_com3;
x_com4 = GEOMETRY.x_com4;
x_s2 = GEOMETRY.x_s2;
x_i = GEOMETRY.x_i;
x_fin = GEOMETRY.x_fin;
x_end=GEOMETRY.x_end;
m=m1+m2+m3+m4;

%% 1.

Wx1=(Wx/m)*m1;
Wz1=(Wz/m)*m1;
Fx1 =(Fx/m)*m1;
Fz1 = (Fz/m)*m1;
L1 = L/3;

P1 = D-Fx1 +Wx1;
R1=  L1 + Wz1 - Fz1;

M1 = ((Fz1+Wz1)*(b1-x_com1)) - L1*(b1-x_com1);
%% 2
Wx2=(Wx/m)*m2;
Wz2=(Wz/m)*m2;
Fx2 =(Fx/m)*m2;
Fz2 = (Fz/m)*m2;

P2 = P1-Fx2+Wx2;
R2= R1 +Wz2 -Fz2;

M2 = ((Fz2+Wz2)*(b2/2)) + M1;
%% 3
Wx3=(Wx/m)*m3;
Wz3=(Wz/m)*m3;
Fx3 =(Fx/m)*m3;
Fz3 = (Fz/m)*m3;

P3 = P2-Fx3+Wx3;
R3= R2 +Wz3 -Fz3;

M3 = ((Fz3+Wz3)*(b3/3)) + M2;
%% 4
Wx4=(Wx/m)*m4;
Wz4=(Wz/m)*m4;
Fx4 =(Fx/m)*m4;
Fz4 =(Fz/m)*m4;
L4 = (2/3)*L;

P4 = P3-Fx4+Wx4;
R4 = R3+Wz4-Fz4;
M4 = ((Fz4+Wz4)*(b4/4)) + M3;

R4_5 = R4+L4;

P5 = P4+Tx;
R5 = R4_5+Tz;
M5 = M4 - L4*(b1+b2+b3+b4-x_fin);


%% PLOT:

% Create the x and y vectors for the plot
x = [0, x_com1,  x_com2, x_com3, x_com, x_com4,x_fin, b1+b2+b3+b4];
y = [D,    P1,     P2,    P3,     P4,  P4,  P4,     P4];

% Plot the piecewise function
figure;
stairs(x, y, 'b-', 'LineWidth', 2); % 'k-' specifies a black line
hold on;

x_points = [ x_com1,  x_com2, x_com3, x_com, x_com4,x_fin, b1+b2+b3+b4];
P_forces = [P1,     P2,    P3,     P4,  P4,  P4,     P4];
for i = 1:length(x_points)
    line([x_points(i), x_points(i)], [D, P_forces(i)], 'Color', 'b', 'LineStyle', '--');
end

% Add labels and markers
xlabel('x');
ylabel('P');
title('Axial Force Distribution');
grid on;
xlim([0, b1+b2+b3+b4]);

%%

% Create the x and y vectors for the plot
x = [0, x_com1,  x_com2, x_com3, x_com, x_com4,x_fin, b1+b2+b3+b4];
y = [L,    R1,     R2,    R3,     R4,  R4,  R4,     R4];

% Plot the piecewise function
figure;
stairs(x, y, 'b-', 'LineWidth', 2); % 'k-' specifies a black line
hold on;

x_points = [ x_com1,  x_com2, x_com3, x_com, x_com4,x_fin, b1+b2+b3+b4];
P_forces = [R1,     R2,    R3,     R4,  R4,  R4,     R4];
for i = 1:length(x_points)
    line([x_points(i), x_points(i)], [L, P_forces(i)], 'Color', 'b', 'LineStyle', '--');
end

% Add labels and markers
xlabel('x');
ylabel('R');
title('Shear Force Distribution');
grid on;
xlim([0, b1+b2+b3+b4]);

%% 

% Create the x and y vectors for the plot
x = [0, x_com1, b1 ,x_s2,x_i ,x_fin, b1+b2+b3+b4];
y = [0,    0,  M1   ,M2,  M3,   M4,  M5];

% Plot the piecewise function
figure;
plot(x, y, 'b-', 'LineWidth', 2); % 'k-' specifies a black line
hold on;

x_points = [0, x_com1, b1 ,x_s2,x_i ,x_fin, b1+b2+b3+b4 ];
P_forces = [0,  M1   ,M2,  M3,   M4,  M5];
for i = 1:length(x_points)
    line([x_points(i), x_points(i)], [L, P_forces(i)], 'Color', 'b', 'LineStyle', '--');
end

% Add labels and markers
xlabel('x');
ylabel('M');
title('Moment Distribution');
grid on;
xlim([0, b1+b2+b3+b4]);

fprintf('Cazzo duro');

end