function [FORCES] = Body_frame_Forces(OUT)

A_11 = cos(OUT.rot_angle);
A_12 = -sin(OUT.rot_angle);
A_21 = sin(OUT.rot_angle);
A_22 = cos(OUT.rot_angle);
D = OUT.D;
gamma = OUT.gamma;
L=OUT.L;
T = OUT.T;
delta = OUT.delta;
theta = OUT.theta;
g = OUT.g;
m = OUT.m;

FORCES.D_x = -A_22*D*sin(pi/2-gamma) - A_21*D*cos(pi/2-gamma);

FORCES.D_z = -A_11*D*cos(pi/2-gamma) - A_12*D*sin(pi/2-gamma);

FORCES.L_x = -A_22*L*cos(pi/2-gamma) + A_21*L*sin(pi/2-gamma);

FORCES.L_z = -A_12*L*cos(pi/2-gamma) + A_11*L*sin(pi/2-gamma);

FORCES.T_x = A_22*T*cos(delta)*cos(theta) + A_21*T*cos(delta)*sin(theta);

FORCES.T_z = A_12*T*cos(delta)*cos(theta) + A_11*T*cos(delta)*cos(theta);

FORCES.G_x = -A_21*g*m;

FORCES.G_z = -A_11*g*m;

FORCES.M = OUT.moment;

FORCES.X_balance =abs(FORCES.D_x + FORCES.L_x + FORCES.T_x + FORCES.G_x);

if FORCES.X_balance>1

    fprintf('Error x balance');


end