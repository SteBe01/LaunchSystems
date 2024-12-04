function [BALANCE]= Load_Distr(OUT,GEOMETRY)

% switch CASE
% 
%     case 2 % Pitch-up 

X_balance = OUT.D + OUT.T_x - OUT.F_x;
Z_balance = OUT.L + OUT.T_z - OUT.F_z + OUT.W;

D=OUT.D;
L=OUT.L;
Tx=OUT.T_x;
Tz=OUT.T_z;
Fx=OUT.F_x;
Fz=OUT.F_z;
W=OUT.W;

if X_balance~=0 && Z_balance~=0
printf('Error in Balance');
end

b1 = GEOMETRY.b1;
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
x_fin = GEOMETRY.x_fin;
m=m1+m2+m3+m4;

%% 1.

W1=(W/m)*m1;
Fx1 =(Fx/m)*m1;
Fz1 = (Fz/m)*m1;
L1 = L/3;

P1 = D-Fx1;
R1=  L1 + W1 - Fz1;
%% 2
W2=(W/m)*m2;
Fx2 =(Fx/m)*m2;
Fz2 = (Fz/m)*m2;

P2 = P1-Fx2;
R2= R1 +W2 -Fz2;
%% 3
W3=(W/m)*m3;
Fx3 =(Fx/m)*m3;
Fz3 = (Fz/m)*m3;

P3 = P2-Fx3;
R3= R2 +W3 -Fz3;
%% 4
W4=(W/m)*m4;
Fx4 =(Fx/m)*m4;
Fz4 =(Fz/m)*m4;
L4 = (2/3)*L;

P4 = P3-Fx4;
R4 = R3+W4-Fz4;

R4_5 = R4+L4;

P5 = P4+Tx;
R5 = R4_5+Tz;




% 
% % Create the x and y vectors for the plot
% x = [0, x_com1, x_com2, x_com3, x_com, x_com4,x_fin,b1+b2+b3+b4];
% y = [D, P1, P2, P3, P4,P4,P4,P5];
% 
% % Plot the piecewise function
% figure;
% stairs(x, y, 'k-', 'LineWidth', 2); % 'k-' specifies a black line
% hold on;
% 
% % Labeling
% xlabel('x');
% ylabel('P');
% title('Piecewise Pressure Function');
% grid on;

% % Annotate x points
% xticks([0, x1, x2, x3, x4]);
% xticklabels({'0', 'x1', 'x2', 'x3', 'x4'});
% ylim([0 max([P1, P2, P3, P4])*1.2]);
% hold off;


% % Step Plot for Shear Forces (R)
% % Segment lengths
% x_segments = [b1, b2, b3, b4];
% x = [0, cumsum(x_segments)]; % Cumulative lengths for plottings
% P_forces = [R_2, R_3, R_4, R_4];
% figure;
% hold on;
% for i = 1:length(x_segments)
%     plot([x(i+1), x(i+1)], [0, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
% end
% stairs(x, [R_1, P_forces], 'LineWidth', 2, 'Color', 'b'); % Step-style plot
% xlabel('x [m]');
% ylabel('Shear Force [N]');
% title('Shear Force Diagram for Pitch-up Maneuver');
% xlim([0, b1+b2+b3+b4]);
% grid on;
% 
% % Segment lengths
% x_segments = [b1, b2, b3, b4]; % Lengths of individual segments
% x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting
% 
% % Moments (M) at segment ends
% M_forces = [M_1, M_2, M_3, M_4]; % Moments at the end of each segment
% 
% % Plot
% figure;
% plot(x, [0, M_forces], 'LineWidth', 2, 'Color', 'b'); % Include zero at start
% xlabel('x [m]');
% ylabel('Moment [Nm]');
% title('Moment Diagram for Pitch-up Maneuver');
% xlim([0, b1+b2+b3+b4]);
% grid on;
% 
% % Add dashed vertical lines for segment boundaries
% hold on;
% for i = 1:length(x_segments)
%     plot([x(i+1), x(i+1)], [0, M_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines


    %case 1 % MAX-Q
% %% 1.
% 
% W1=(W/m)*m1;
% Fx1 =(Fx/m)*m1;
% P1 = D-Fx1;
% %% 2
% W2=(W/m)*m2;
% Fx2 =(Fx/m)*m2;
% P2 = P1-Fx2;
% %% 3
% W3=(W/m)*m3;
% Fx3 =(Fx/m)*m3;
% P3 = P2-Fx3;
% %% 4
% W4=(W/m)*m4;
% Fx4 =(Fx/m)*m4;
% P4 = P3-Fx4;
% P5 = P3-Fx4+Tx;

%end % switch
end % function