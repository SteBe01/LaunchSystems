function [GEOMETRY] = GEO(h_it,M_it)

GEOMETRY.M_pay = M_it.pay;
GEOMETRY.OF1 = M_it.stg1.OF;
GEOMETRY.OF2 = M_it.stg1.OF;
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

M_w_d =  M_it.wiring/10;
M_a_d =  M_it.avionics/10;

GEOMETRY.M_s_1 = M_it.str1;
GEOMETRY.M_s_2 = M_it.str2;

GEOMETRY.m1=M_it.pay_effective+M_it.stg1.fairing;      % fairing
GEOMETRY.m2=M_w_d+M_a_d+M_it.stg2.C1.m;                  % forward skirt
GEOMETRY.m3=M_it.stg2.tot_rp1+M_w_d+M_a_d;             % fuel tank2
%GEOMETRY.m34 = M_it.stg1.C2.m+M_w_d+M_a_d;               % intertank 2
GEOMETRY.m4=M_it.stg2.tot_lox+M_w_d+M_a_d+M_it.stg2.lox_insulation;  % oxygen tank2
GEOMETRY.m5=M_it.stg2.motor+M_w_d+M_a_d+M_it.stg2.T_struct+M_it.stg2.C3.m; % aft skirt 2
GEOMETRY.m6=+M_w_d+M_a_d+M_it.stg1.recovery+M_it.stg1.C1.m;            % interstage
GEOMETRY.m7=M_it.stg1.tot_rp1+M_w_d+M_a_d;               %  fuel tank1
%GEOMETRY.m78 = M_it.stg1.C2.m+M_a_d+M_w_d;                % intertank 1
GEOMETRY.m8=M_it.stg1.tot_lox+M_w_d+M_a_d+M_it.stg1.lox_insulation;   %  oxygen tank1
GEOMETRY.m9=M_w_d+M_a_d+M_it.stg1.T_struct+M_it.stg1.C3.m;           % aft skirt 1
GEOMETRY.m10=M_it.stg1.motor+M_w_d+M_a_d;               % engine

DELTA = (GEOMETRY.M01-GEOMETRY.m1-GEOMETRY.m2-GEOMETRY.m3-GEOMETRY.m4-GEOMETRY.m5-GEOMETRY.m6-GEOMETRY.m7-GEOMETRY.m8-GEOMETRY.m9-GEOMETRY.m10);

DELTA1 = 0.5*DELTA;
DELTA2 = 0.2*DELTA;
DELTA3 = 0.3*DELTA;

GEOMETRY.m1=M_it.pay_effective+M_it.stg1.fairing;      % fairing
GEOMETRY.m2=M_w_d+M_a_d+M_it.stg2.C1.m+DELTA1;                  % forward skirt
GEOMETRY.m3=M_it.stg2.tot_rp1+M_w_d+M_a_d;             % fuel tank2
%GEOMETRY.m34 = M_it.stg1.C2.m+M_w_d+M_a_d+DELTA;               % intertank 2
GEOMETRY.m4=M_it.stg2.tot_lox+M_w_d+M_a_d+M_it.stg2.lox_insulation;  % oxygen tank2
GEOMETRY.m5=M_it.stg2.motor+M_w_d+M_a_d+M_it.stg2.T_struct+M_it.stg2.C3.m+DELTA2; % aft skirt 2
GEOMETRY.m6=+M_w_d+M_a_d+M_it.stg1.recovery+M_it.stg1.C1.m;            % interstage
GEOMETRY.m7=M_it.stg1.tot_rp1+M_w_d+M_a_d;                   %  fuel tank1
%GEOMETRY.m78 = M_it.stg1.C2.m+M_a_d+M_w_d+DELTA;                  % intertank 1
GEOMETRY.m8=M_it.stg1.tot_lox+M_w_d+M_a_d+M_it.stg1.lox_insulation;   %  oxygen tank1
GEOMETRY.m9=M_w_d+M_a_d+M_it.stg1.T_struct+M_it.stg1.C3.m+DELTA3;           % aft skirt 1
GEOMETRY.m10=M_it.stg1.motor+M_w_d+M_a_d;               % engine

GEOMETRY.m_lox_2 = M_it.stg2.prop * M_it.stg2.OF / (1+M_it.stg2.OF);%[kg] mass of lox
GEOMETRY.m_fuel_2 =  M_it.stg2.prop * 1  / (1+M_it.stg2.OF);%[kg] mass of rp1
GEOMETRY.m_prop_2 =  M_it.stg2.prop;
GEOMETRY.m_tank_f_2 = M_it.stg2.tank_rp1;
GEOMETRY.m_tank_lox_2 = M_it.stg2.tank_lox;

GEOMETRY.m_lox_1 = M_it.stg1.prop * M_it.stg1.OF / (1+M_it.stg1.OF);%[kg] mass of lox
GEOMETRY.m_fuel_1 =  M_it.stg1.prop * 1  / (1+M_it.stg1.OF);%[kg] mass of rp1
GEOMETRY.m_prop_1 =  M_it.stg1.prop;
GEOMETRY.m_tank_f_1 = M_it.stg1.tank_rp1;
GEOMETRY.m_tank_lox_1 = M_it.stg1.tank_lox;

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

GEOMETRY.eps_1 =abs((GEOMETRY.m6+GEOMETRY.m7+GEOMETRY.m8+GEOMETRY.m9+GEOMETRY.m10-GEOMETRY.m_prop_1)) /(GEOMETRY.m6+GEOMETRY.m7+GEOMETRY.m8+GEOMETRY.m9+GEOMETRY.m10);

GEOMETRY.eps_2 =abs((GEOMETRY.m1+GEOMETRY.m2+GEOMETRY.m3+GEOMETRY.m4+GEOMETRY.m5-GEOMETRY.m_prop_2-GEOMETRY.M_pay)) /(GEOMETRY.m1+GEOMETRY.m2+GEOMETRY.m3+GEOMETRY.m4+GEOMETRY.m5-GEOMETRY.M_pay);


end