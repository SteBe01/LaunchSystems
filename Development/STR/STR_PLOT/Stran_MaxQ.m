function [THICKNESS]= Stran_MaxQ(FORCES,GEOMETRY,MAT,OUT)

D = FORCES.D_x + FORCES.L_x;
L = FORCES.D_z + FORCES.L_z;

T_x = FORCES.T_x;
T_z = FORCES.T_z;

F_x = -FORCES.G_x + FORCES.F_x;
F_z = -FORCES.G_z + FORCES.F_z;

FS = OUT.FS;

sigma_yield_2 = MAT.MAT_2.sy;
sigma_yield_1 = MAT.MAT_1.sy;

b1 = GEOMETRY.Stk1.Cone.L;

b2 = GEOMETRY.Stk1.Stage2.L;

b3 = GEOMETRY.Stk1.Interstage.L;

b4 = GEOMETRY.Stk1.Stage1.L;

% m1 =  GEOMETRY.Stk1.Cone.m_cone;
% 
% m2 = GEOMETRY.Stk1.Stage2.m_stage2;
% 
% m3 = GEOMETRY.Stk1.Interstage.m_interstage;
% 
% m4 = GEOMETRY.Stk1.Stage1.m_stage1;

m1 = 190*OUT.m/10000;
m2 = 1089*OUT.m/10000;
m3 = 181*OUT.m/10000;
m4 = 8540*OUT.m/10000;

Diam_1 = GEOMETRY.Diam_1;
Diam_2 = GEOMETRY.Diam_2;

l_tot = GEOMETRY.Length;

%% Load division:

q_D = D/l_tot;
q_L = L/l_tot;

n_x = F_x/(m1+m2+m3+m4);
n_z = F_z/(m1+m2+m3+m4);




end