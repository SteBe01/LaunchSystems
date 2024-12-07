function [GEOMETRY] = GEO(h_it,M_it,M1,M2)


GEOMETRY.M_pay = M_it.pay;

% GEOMETRY.l_tot= h_it.tot;
% GEOMETRY.b1 = h_it.fairing;  % fairing length
% GEOMETRY.b2 = h_it.stg2.C1;  % forward skirt length
% GEOMETRY.b3 = h_it.stg2.cyl_rp1;  % height of fuel tank2
% GEOMETRY.b34 = h_it.stg2.C2; % 16 cm
% GEOMETRY.b4 = h_it.stg2.cyl_lox; %  height of oxygen tank2
% GEOMETRY.b5 = h_it.stg2.C3;%  aft skirt length 2
% GEOMETRY.b6 = h_it.stg1.C1; % interstage length
% GEOMETRY.b7 = h_it.stg1.cyl_rp1;% height of fuel tank1
% GEOMETRY.b78 =h_it.stg1.C2; % 16 cm
% GEOMETRY.b8 = h_it.stg1.cyl_lox;% height of oxygen tank1
% GEOMETRY.b9 = h_it.stg1.C3; % aft skirt length 1
% GEOMETRY.b10 = h_it.stg1.motor; % h engine
GEOMETRY.l_tot= h_it.tot;
GEOMETRY.l_tot_2 = h_it.stg2.tot;
GEOMETRY.b1 = h_it.fairing;  % fairing length
GEOMETRY.b2 = h_it.stg2.C1;  % forward skirt length
GEOMETRY.b3 = h_it.stg2.cyl_rp1;  % height of fuel tank2
GEOMETRY.b34 = h_it.stg2.C2; % 16 cm
GEOMETRY.b4 = h_it.stg2.cyl_lox; %  height of oxygen tank2
GEOMETRY.b5 = h_it.stg2.C3;%  aft skirt length 2
GEOMETRY.b6 = h_it.stg1.C1; % interstage length
GEOMETRY.b7 = h_it.stg1.cyl_rp1;% height of fuel tank1
GEOMETRY.b78 =h_it.stg1.C2; % 16 cm
GEOMETRY.b8 = h_it.stg1.cyl_lox;% height of oxygen tank1
GEOMETRY.b9 = h_it.stg1.C3; % aft skirt length 1
GEOMETRY.b10 = h_it.stg1.motor; % h engine

GEOMETRY.Diam_1 = M_it.diam1;
GEOMETRY.Diam_2 = M_it.diam2;

GEOMETRY.M01 = M_it.M0;

GEOMETRY.x_com = GEOMETRY.l_tot - h_it.CG;
GEOMETRY.X_COM2=GEOMETRY.l_tot-h_it.stg2.CG.tot;

GEOMETRY.M_avionics = M_it.avionics;
GEOMETRY.M_cables = M_it.wiring;

GEOMETRY.M_s_1 = M_it.str1;
GEOMETRY.M_s_2 = M_it.str2;

M_remaining_1 = M1.tot-M1.tot_rp1-M1.tot_lox-M1.motor;
M_remaining_2 = M2.tot-M_it.adapter-M2.tot_rp1-M2.tot_lox-M2.motor;
M_distr=M_it.wiring+ M_it.avionics/11;
M_diff =  97.9259/10;
Mdiff_2 = -150.9240/10;
GEOMETRY.m1=M_it.pay_effective+M_distr+Mdiff_2;  % fairing
GEOMETRY.m2=M_remaining_2+M_distr-M_diff+Mdiff_2;% forward skirt
GEOMETRY.m3=M2.tot_rp1-M_distr-M_diff+Mdiff_2; % fuel tank2
GEOMETRY.m4=M2.tot_lox-M_distr-M_diff+Mdiff_2;% oxygen tank2
GEOMETRY.m5=M2.motor+M_distr-M_diff+Mdiff_2; %  aft skirt 2
GEOMETRY.m6=M_remaining_1/2 -M_distr-M_diff+Mdiff_2;% interstage
GEOMETRY.m7=M1.tot_rp1+M_distr-M_diff-9.4577+Mdiff_2;%  fuel tank1
GEOMETRY.m8=M1.tot_lox+M_distr-M_diff+Mdiff_2;%  oxygen tank1
GEOMETRY.m9=M_remaining_1/2 +M_distr-M_diff+Mdiff_2;% aft skirt 1
GEOMETRY.m10=M1.motor+M_distr-M_diff+Mdiff_2; % engine

% GEOMETRY.m_lox_2 = M2.prop * M2.OF / (1+M2.OF);%[kg] mass of lox
% GEOMETRY.m_fuel_2 =  M2.prop * 1  / (1+M2.OF);%[kg] mass of rp1
% GEOMETRY.m_prop_2 = M2.prop;
GEOMETRY.m_tank_f_2 = M2.tank_rp1;
GEOMETRY.m_tank_lox_2 = M2.tank_lox;

% GEOMETRY.m_lox_1 = M1.prop * M1.OF / (1+M1.OF);%[kg] mass of lox
% GEOMETRY.m_fuel_1 =  M1.prop * 1  / (1+M1.OF);%[kg] mass of rp1
% GEOMETRY.m_prop_1 = M1.prop;
GEOMETRY.m_tank_f_1 = M1.tank_rp1;
GEOMETRY.m_tank_lox_1 = M1.tank_lox;

STAGE1=1;
[J_yaw1,J_pitch1,J_roll1,J_yaw_empty1,J_pitch_empty1,J_roll_empty1]= MOI(GEOMETRY,STAGE1);
STAGE2=2;
[J_yaw2,J_pitch2,J_roll2,J_yaw_empty2,J_pitch_empty2,J_roll_empty2]= MOI(GEOMETRY,STAGE2);

GEOMETRY.Jyaw1 = J_yaw1;
GEOMETRY.Jyaw2 = J_yaw2;
GEOMETRY.Jpitch1 = J_pitch1;
GEOMETRY.Jpitch2 = J_pitch2;
GEOMETRY.Jroll2 = J_roll2;
GEOMETRY.Jroll1 = J_roll1;

GEOMETRY.Jyaw_empty1 = J_yaw_empty1;
GEOMETRY.Jyaw_empty2= J_yaw_empty2;
GEOMETRY.Jpitch_empty_1=J_pitch_empty1;
GEOMETRY.Jpitch_empty_2=J_pitch_empty2;
GEOMETRY.Jroll_empty1 = J_roll_empty1;
GEOMETRY.Jroll_empty2 = J_roll_empty2;

GEOMETRY.thick_1_max = max([M_it.th1.C1,M_it.th1.C2,M_it.th1.C3,M_it.th1.rp1,M_it.th1.lox]);
GEOMETRY.thick_2_max = max([M_it.th2.C1,M_it.th2.C2,M_it.th2.C3,M_it.th2.rp1,M_it.th2.lox]);
end