function [GEOMETRY] = GEO(h_it,M_it,M1,M2)


GEOMETRY.M_pay = M_it.pay;

GEOMETRY.l_tot= h_it.tot;
GEOMETRY.b1 = h_it.fairing;  % fairing length
GEOMETRY.b2 = h_it.stg2.C1;  % forward skirt length
GEOMETRY.b3 = h_it.stg2.tank_rp1;  % height of fuel tank2
GEOMETRY.b34 = h_it.stg2.C2; % 16 cm
GEOMETRY.b4 = h_it.stg2.tank_lox; %  height of oxygen tank2
GEOMETRY.b5 = h_it.stg2.til_tank;%  aft skirt length 2
GEOMETRY.b6 = h_it.stg1.C1; % interstage length
GEOMETRY.b7 = h_it.stg1.tank_rp1;% height of fuel tank1
GEOMETRY.b78 =h_it.stg1.C2; % 16 cm
GEOMETRY.b8 = h_it.stg1.tank_lox;% height of oxygen tank1
GEOMETRY.b9 = h_it.stg1.C3; % aft skirt length 1
GEOMETRY.b10 = h_it.stg1.motor; % h engine

GEOMETRY.Diam_1 = M_it.diam1;
GEOMETRY.Diam_2 = M_it.diam2;

GEOMETRY.M01 = M_it.M0;

GEOMETRY.x_com = h_it.CG;

GEOMETRY.M_avionics = M_it.avionics;
GEOMETRY.M_cables = M_it.wiring;

GEOMETRY.M_s_1 = M_it.str1;
GEOMETRY.M_s_2 = M_it.str2;

M_remaining = M_it.M0 - M_it.adapter -M2.tot_rp1-M2.tot_lox- M2.motor-M1.tot_rp1-M1.tot_lox-M1.motor;

GEOMETRY.m1=M_it.adapter;  % fairing
%GEOMETRY.m2=;% forward skirt
 GEOMETRY.m3=M2.tot_rp1; % fuel tank2
GEOMETRY.m4=M2.tot_lox;% oxygen tank2
GEOMETRY.m5=M2.motor; %  aft skirt 2
% GEOMETRY.m6=;% interstage
GEOMETRY.m7=M1.tot_rp1;%  fuel tank1
GEOMETRY.m8=M1.tot_lox;%  oxygen tank1
%GEOMETRY.m9=;% aft skirt 1
GEOMETRY.m10=M1.motor; % engine




end