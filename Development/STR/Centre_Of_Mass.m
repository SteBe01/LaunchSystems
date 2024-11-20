function [COM,MoI] = Centre_Of_Mass(TR,Tank1,Tank2,Engine1,Engine2,t)

%% Main outputs:

% COM.TR.x_LV = C.O.M. position wrt nose (x=0) at t = t1+t2
% COM.TR.x_LV0 = C.O.M. position wrt nose (x=0) at t = 0

% MoI.TR.Jy = Inertia pitch moment at  t = t1+t2
% MoI.TR.Jy0 = Inertia pitch moment at  t = 0

% structures to be defined for Engine and t
% t.t1 = 0;
% t.t2 = 0;
% 
% t.t_burn1 = 700;
% t.t_burn2 = 500;
%
% Engine2.L_nozzle = L_nozzle2;
% Engine1.L_nozzle= L_nozzle1;
% 
% Engine2.m_dot_lox = m_dotlox_2;
%  Engine2.m_dot_fuel=m_dotfuel_2;
% Engine1.m_dot_lox=m_dotlox_1;
% Engine1.m_dot_fuel=m_dotfuel_1;
% 
% Engine2.M = 0;
% Engine1.M = 0;
% Engine2.t_2 = 0;
% Engine1.t_1 = 0;
% Engine1.R_ext = 1;
% Engine2.R_ext = 1;

L_nozzle2 = Engine2.L_nozzle;
L_nozzle1 = Engine1.L_nozzle;

m_dotlox_2 = Engine2.m_dot_lox;
m_dotfuel_2 = Engine2.m_dot_fuel;
m_dotlox_1 = Engine1.m_dot_lox;
m_dotfuel_1 = Engine1.m_dot_fuel;

if strcmp(Tank2.OX.FLAG, 'Cyl') && strcmp(Tank2.FU.FLAG, 'Cyl')
Diam_2 = TR.Diameter;

elseif strcmp(Tank2.OX.FLAG, 'Sphere') && strcmp(Tank2.FU.FLAG, 'Sphere')
Diam_2 = 2*(max([Tank2.OX.R_sphere_lox;Tank2.FU.R_sphere_fuel]));

end

if strcmp(Tank1.OX.FLAG, 'Cyl') && strcmp(Tank1.FU.FLAG, 'Cyl')
Diam_1 = TR.Diameter;

elseif strcmp(Tank1.OX.FLAG, 'Sphere') && strcmp(Tank1.FU.FLAG, 'Sphere')
Diam_1 = 2*(max([Tank1.OX.R_sphere_lox;Tank1.FU.R_sphere_fuel]));

end

M_cables = 1.058*(TR.Length^(0.25))*sqrt(TR.M.M01); % [kg] Empirical formula slides Maggi 06, structures part 1

M_avionics = 10*(TR.M.M01)^(0.361); % [kg] Empirical formula slides Maggi 06, structures part 1

M_cables_distributed = M_cables/19;

% To distinguish between 1st part of flight and 2nd part after stage 1 is
% jettisoned

if t.t1 <t.t_burn1

%% STACK 1: t = [0-t1]

Ms1_stage_ratio = 1- TR.M.Ms2/TR.M.Ms1;

M_avionics_s1 = M_avionics*Ms1_stage_ratio;

Ms2_stage_ratio = TR.M.Ms2/TR.M.Ms1;

M_avionics_s2 = M_avionics*Ms2_stage_ratio;

t1 = t.t1;
t2 = 0; % tandem launcher, so no prop of 2nd satge is burned

% x = 0 at the tip

% First of all of the lengths are defined as in slide 10 of str mass sildes
% Percentage of 5% added to L fair and 10% to interstage to account for
% space to place parachute

% 2 stage

x_fair = TR.Fair.fair_length*(1+0.05);

if strcmp(Tank2.OX.FLAG, 'Cyl') && strcmp(Tank2.FU.FLAG, 'Cyl')
% Top FU, bottom OX
L_fwrd_skirt_s2 = ((1/3)*Diam_2) + Tank2.FU.H_dome;
L_cyl_tank_FU2 = Tank2.FU.H_cyl;
L_connection_2 = ((1/4)*Diam_2) + (Tank2.OX.H_dome + Tank2.FU.H_dome);
L_cyl_tank_OX2 = Tank2.OX.H_cyl;
L_dwrd_skirt_s2 = ((1/3)*Diam_2) + Tank2.OX.H_dome;

x_fwrd_skirt_s2 = x_fair + L_fwrd_skirt_s2;
x_cyl_22 = x_fwrd_skirt_s2 + L_cyl_tank_FU2;
x_connection_2 = x_cyl_22 + L_connection_2;
x_cyl_21 = x_connection_2 + L_cyl_tank_OX2;
x_dwrd_skirt_s2 = x_cyl_21 + L_dwrd_skirt_s2;

elseif strcmp(Tank2.OX.FLAG, 'Sphere') && strcmp(Tank2.FU.FLAG, 'Sphere')

   L_fwrd_skirt_s2 = (Tank2.FU.R_sphere_fuel/2)*(1+0.05);
   L_rem_sphere_FU2 = Tank2.FU.R_sphere_fuel;
   L_connection_2 = ((1/4)*Diam_2) + ((Tank2.OX.R_sphere_lox/2) + (Tank2.FU.R_sphere_fuel/2));
   L_rem_sphere_OX2 = Tank2.OX.R_sphere_lox;
   L_dwrd_skirt_s2 = (Tank2.OX.R_sphere_lox/2) + ((1/3)*Diam_2);

x_fwrd_skirt_s2 = x_fair + L_fwrd_skirt_s2;
x_sph_22 = x_fwrd_skirt_s2 + L_rem_sphere_FU2;
x_connection_2 = x_sph_22 + L_connection_2;
x_sph_21 = x_connection_2 + L_rem_sphere_OX2;
x_dwrd_skirt_s2 = x_sph_21 + L_dwrd_skirt_s2;


end

L_interstage = (5/4)*Diam_1*(1+0.1);%L_nozzle2*(1+0.1);

x_interstage = x_dwrd_skirt_s2 + L_interstage;

% 1 stage

if strcmp(Tank1.OX.FLAG, 'Cyl') && strcmp(Tank1.FU.FLAG, 'Cyl')
% Top FU, bottom OX
L_cyl_tank_FU1 = Tank1.FU.H_cyl;
L_connection_1 = ((1/4)*Diam_1) + (Tank1.OX.H_dome + Tank1.FU.H_dome);
L_cyl_tank_OX1 = Tank1.OX.H_cyl;
L_end = Diam_1 + L_nozzle1;

x_cyl_12 = x_interstage + L_cyl_tank_FU1;
x_connection_1 = x_cyl_12 + L_connection_1;
x_cyl_11 = L_cyl_tank_OX1 + x_connection_1;
x_end = x_cyl_11 + L_end;

elseif strcmp(Tank1.OX.FLAG, 'Sphere') && strcmp(Tank1.FU.FLAG, 'Sphere')

   L_rem_sphere_FU1 = (Tank1.FU.R_sphere_fuel);
   L_connection_1 = ((1/4)*Diam_1) + ((Tank1.OX.R_sphere_lox/2) + (Tank1.FU.R_sphere_fuel/2));
   L_rem_sphere_OX1 = Tank1.OX.R_sphere_lox;
   L_end = Diam_1 + L_nozzle1; 

x_sph_12 = x_interstage + L_rem_sphere_FU1;
x_connection_1 = x_sph_12 + L_connection_1;
x_sph_11 = L_rem_sphere_OX1 + x_connection_1;
x_end = x_sph_11 + L_end;

end

Delta_length_err = x_end - TR.Length;

if x_end > TR.Length

  fprintf('Length overbudget: \n Delta_length_err = %d [m] \n ',abs(Delta_length_err));

x_end = TR.Length;

if strcmp(Tank1.OX.FLAG, 'Cyl') && strcmp(Tank1.FU.FLAG, 'Cyl') % correct end length

    L_end = TR.Length - x_cyl_11;

elseif strcmp(Tank1.OX.FLAG, 'Sphere') && strcmp(Tank1.FU.FLAG, 'Sphere')
    
L_end = TR.Length - x_sph_11;

end

end


%% C.O.M. positioning for each component

% x=0 at tip

COM.Stk1.Pay.x_pay = TR.Fair.fair_length/2;

COM.Stk1.Fair.x_fair = TR.Fair.fair_length*(2/3);

COM.Stk1.Prac2.x_parachute_2 = x_fair + L_fwrd_skirt_s2/2;

if strcmp(Tank2.OX.FLAG, 'Cyl') && strcmp(Tank2.FU.FLAG, 'Cyl')

COM.Stk1.Prop2_FU.x_Mp2_FU0 = x_fwrd_skirt_s2+ L_cyl_tank_FU2/2;
COM.Stk1.Prop2_FU.x_Mp2_FU = COM.Stk1.Prop2_FU.x_Mp2_FU0 - (((m_dotfuel_2*t2)/(Tank2.FU.M_fuel))*COM.Stk1.Prop2_FU.x_Mp2_FU0);
COM.Stk1.Tank2_FU.x_tank2_FU = COM.Stk1.Prop2_FU.x_Mp2_FU0;

COM.Stk1.Prop2_OX.x_Mp2_OX0 = x_connection_2 + L_cyl_tank_OX2/2;
COM.Stk1.Prop2_OX.x_Mp2_OX = COM.Stk1.Prop2_OX.x_Mp2_OX0 - (((m_dotlox_2*t2)/(Tank2.OX.M_lox))*COM.Stk1.Prop2_OX.x_Mp2_OX0);
COM.Stk1.Tank2_OX.x_tank2_OX = COM.Stk1.Prop2_OX.x_Mp2_OX0;

COM.Stk1.Pipes2.x_pipes2 = x_cyl_22 + Tank2.FU.H_dome; % about interface between 2 tanks

COM.Stk1.Aft_Struct2.x_aft_struct2 = x_cyl_21 + L_dwrd_skirt_s2/2;

COM.Stk1.Engine2.x_engine2 = x_cyl_21 + ((L_dwrd_skirt_s2)*(2/3));

elseif strcmp(Tank2.OX.FLAG, 'Sphere') && strcmp(Tank2.FU.FLAG, 'Sphere')

COM.Stk1.Prop2_FU.x_Mp2_FU0 =  x_fwrd_skirt_s2+ L_rem_sphere_FU2/2;
COM.Stk1.Prop2_FU.x_Mp2_FU = COM.Stk1.Prop2_FU.x_Mp2_FU0 - (((m_dotfuel_2*t2)/(Tank2.FU.M_fuel))*COM.Stk1.Prop2_FU.x_Mp2_FU0);
COM.Stk1.Tank2_FU.x_tank2_FU = COM.Stk1.Prop2_FU.x_Mp2_FU0;

COM.Stk1.Prop2_OX.x_Mp2_OX0 = x_connection_2 + L_rem_sphere_OX2/2;
COM.Stk1.Prop2_OX.x_Mp2_OX = COM.Stk1.Prop2_OX.x_Mp2_OX0 - (((m_dotlox_2*t2)/(Tank2.OX.M_lox))*COM.Stk1.Prop2_OX.x_Mp2_OX0);
COM.Stk1.Tank2_OX.x_tank2_OX = COM.Stk1.Prop2_OX.x_Mp2_OX0;

COM.Stk1.Pipes2.x_pipes2 = x_sph_22 + Tank2.FU.R_sphere_fuel/2; % about interface between 2 tanks

COM.Stk1.Aft_Struct2.x_aft_struct2 = x_sph_21 + L_dwrd_skirt_s2/2;

COM.Stk1.Engine2.x_engine2 = x_sph_21 + ((L_dwrd_skirt_s2)*(2/3));

end

COM.Stk1.Interstage.x_interstage = x_dwrd_skirt_s2 + (L_interstage*(1/2)); % Position that accounts for mass of nozzle + parachute + part of avionics; 1/2: to be changed maybe

if strcmp(Tank1.OX.FLAG, 'Cyl') && strcmp(Tank1.FU.FLAG, 'Cyl')

COM.Stk1.Prop1_FU.x_Mp1_FU0 = x_interstage + Tank1.FU.H_cyl/2;
COM.Stk1.Prop1_FU.x_Mp1_FU = COM.Stk1.Prop1_FU.x_Mp1_FU0 - (((m_dotfuel_1*t1)/(Tank1.FU.M_fuel))*COM.Stk1.Prop1_FU.x_Mp1_FU0);
COM.Stk1.Tank1_FU.x_tank1_FU = COM.Stk1.Prop1_FU.x_Mp1_FU0;

COM.Stk1.Pipes1.x_pipes1 = x_cyl_12 + Tank1.FU.H_dome;

COM.Stk1.Prop1_OX.x_Mp1_OX0 = x_connection_1 + Tank1.OX.H_cyl/2;
COM.Stk1.Prop1_OX.x_Mp1_OX = COM.Stk1.Prop1_OX.x_Mp1_OX0 - (((m_dotlox_1*t1)/(Tank1.OX.M_lox))*COM.Stk1.Prop1_OX.x_Mp1_OX0);
COM.Stk1.Tank1_OX.x_tank1_OX=COM.Stk1.Prop1_OX.x_Mp1_OX0;

elseif strcmp(Tank1.OX.FLAG, 'Sphere') && strcmp(Tank1.FU.FLAG, 'Sphere')

COM.Stk1.Prop1_FU.x_Mp1_FU0 = x_interstage + L_rem_sphere_FU1/2;
COM.Stk1.Prop1_FU.x_Mp1_FU = COM.Stk1.Prop1_FU.x_Mp1_FU0 - (((m_dotfuel_1*t1)/(Tank1.FU.M_fuel))*COM.Stk1.Prop1_FU.x_Mp1_FU0);
COM.Stk1.Tank1_FU.x_tank1_FU = COM.Stk1.Prop1_FU.x_Mp1_FU0;


COM.Stk1.Pipes1.x_pipes1 = x_sph_12 + Tank1.FU.R_sphere_fuel/2;

COM.Stk1.Prop1_OX.x_Mp1_OX0 = x_connection_1 + L_rem_sphere_OX1/2;
COM.Stk1.Prop1_OX.x_Mp1_OX = COM.Stk1.Prop1_OX.x_Mp1_OX0 - (((m_dotlox_1*t1)/(Tank1.OX.M_lox))*COM.Stk1.Prop1_OX.x_Mp1_OX0);
COM.Stk1.Tank1_OX.x_tank1_OX=COM.Stk1.Prop1_OX.x_Mp1_OX0;

end

COM.Stk1.Aft_Struct1.x_aft_struct1 = x_cyl_11 + (L_end*(1/4)); % Sustaining structures

COM.Stk1.Engine1.x_engine1 = x_cyl_11 + (L_end*(1/2)); % Engine + nozzle

COM.Stk1.Fins.x_fins = x_cyl_11 + (L_end*(3/4)); % x coord of fins

%% Mass association to every quantity in COM struct:

COM.Stk1.Pay.m_pay = 250;

COM.Stk1.Fair.m_fair = TR.Fair.M_fair +M_cables_distributed;

COM.Stk1.Prac2.m_parachute_2 = 0.1*TR.M.Ms2+ M_cables_distributed + M_avionics_s2/2;

COM.Stk1.Tank2_FU.m_tank2_FU = Tank2.FU.M_tot_tank_fuel+ M_cables_distributed;

COM.Stk1.Tank2_OX.m_tank2_OX = Tank2.OX.M_tot_tank_lox+ M_cables_distributed;

COM.Stk1.Prop2_FU.m_prop2_FU0 = Tank2.FU.M_fuel;
COM.Stk1.Prop2_FU.m_prop2_FU = COM.Stk1.Prop2_FU.m_prop2_FU0 - (m_dotfuel_2*t2);

COM.Stk1.Prop2_OX.m_prop2_OX0 = Tank2.OX.M_lox;
COM.Stk1.Prop2_OX.m_prop2_OX = COM.Stk1.Prop2_OX.m_prop2_OX0 - (m_dotlox_2*t2);

COM.Stk1.Pipes2.m_pipes2 = 10+ M_cables_distributed;

COM.Stk1.Aft_Struct2.m_aft_struct2 = 100+ M_cables_distributed + M_avionics_s2/2;

COM.Stk1.Engine2.m_engine2 = 300+ M_cables_distributed;

COM.Stk1.Interstage.m_interstage = 0.1*TR.M.Ms1 + 20+ M_cables_distributed + M_avionics_s1/2;

COM.Stk1.Tank1_FU.m_tank1_FU = Tank1.FU.M_tot_tank_fuel+ M_cables_distributed;

COM.Stk1.Tank1_OX.m_tank1_OX = Tank1.OX.M_tot_tank_lox+ M_cables_distributed;

COM.Stk1.Prop1_FU.m_prop1_FU0 = Tank1.FU.M_fuel;
COM.Stk1.Prop1_FU.m_prop1_FU = COM.Stk1.Prop1_FU.m_prop1_FU0 - ((m_dotfuel_1*t1));

COM.Stk1.Prop1_OX.m_prop1_OX0 = Tank1.OX.M_lox;
COM.Stk1.Prop1_OX.m_prop1_OX = COM.Stk1.Prop1_OX.m_prop1_OX0 - ((m_dotlox_1*t1));

COM.Stk1.Pipes1.m_pipes1 = 10+ M_cables_distributed;

COM.Stk1.Aft_Struct1.m_aft_struct1 = 300+ M_cables_distributed +M_avionics_s1;

COM.Stk1.Engine1.m_engine1 = 200+ M_cables_distributed; % Engine + nozzle

COM.Stk1.Fins.m_fins = 10+ M_cables_distributed; % x coord of fins

%% Computation of position of COM wrt nose and M0

% Vector assembly ( careful with associating right components)

COM_vec_m = [COM.Stk1.Pay.m_pay;COM.Stk1.Fair.m_fair;COM.Stk1.Prac2.m_parachute_2;COM.Stk1.Prop2_FU.m_prop2_FU;COM.Stk1.Tank2_FU.m_tank2_FU;COM.Stk1.Prop2_OX.m_prop2_OX;COM.Stk1.Tank2_OX.m_tank2_OX;COM.Stk1.Pipes2.m_pipes2;COM.Stk1.Aft_Struct2.m_aft_struct2;COM.Stk1.Engine2.m_engine2;COM.Stk1.Interstage.m_interstage;COM.Stk1.Prop1_FU.m_prop1_FU;COM.Stk1.Tank1_FU.m_tank1_FU;COM.Stk1.Pipes1.m_pipes1;COM.Stk1.Prop1_OX.m_prop1_OX;COM.Stk1.Tank1_OX.m_tank1_OX;COM.Stk1.Aft_Struct1.m_aft_struct1;COM.Stk1.Engine1.m_engine1;COM.Stk1.Fins.m_fins];
COM_vec_m0 = [COM.Stk1.Pay.m_pay;COM.Stk1.Fair.m_fair;COM.Stk1.Prac2.m_parachute_2;COM.Stk1.Prop2_FU.m_prop2_FU0;COM.Stk1.Tank2_FU.m_tank2_FU;COM.Stk1.Prop2_OX.m_prop2_OX0;COM.Stk1.Tank2_OX.m_tank2_OX;COM.Stk1.Pipes2.m_pipes2;COM.Stk1.Aft_Struct2.m_aft_struct2;COM.Stk1.Engine2.m_engine2;COM.Stk1.Interstage.m_interstage;COM.Stk1.Prop1_FU.m_prop1_FU0;COM.Stk1.Tank1_FU.m_tank1_FU;COM.Stk1.Pipes1.m_pipes1;COM.Stk1.Prop1_OX.m_prop1_OX0;COM.Stk1.Tank1_OX.m_tank1_OX;COM.Stk1.Aft_Struct1.m_aft_struct1;COM.Stk1.Engine1.m_engine1;COM.Stk1.Fins.m_fins];
%COM_vec_m0 = [COM.Pay.m_pay;COM.Fair.m_fair;COM.Prac2.m_parachute_2;COM.Prop2_FU.m_prop2_FU0;COM.Tank2_FU.m_tank2_FU;COM.Prop2_OX.m_prop2_OX0;COM.Tank2_OX.m_tank2_OX;COM.Pipes2.m_pipes2;COM.Aft_Struct2.m_aft_struct2;COM.Engine2.m_engine2;COM.Interstage.m_interstage;COM.Prop1_FU.m_prop1_FU0;COM.Tank1_FU.m_tank1_FU;COM.Pipes1.m_pipes1;COM.Prop1_OX.m_prop1_OX0;COM.Tank1_OX.m_tank1_OX;COM.Stk1.Aft_Struct1.m_aft_struct1;COM.Stk1.Engine1.m_engine1;COM.Stk1.Fins.m_fins];
COM_vec_x = [COM.Stk1.Pay.x_pay;COM.Stk1.Fair.x_fair;COM.Stk1.Prac2.x_parachute_2;COM.Stk1.Prop2_FU.x_Mp2_FU;COM.Stk1.Tank2_FU.x_tank2_FU;COM.Stk1.Prop2_OX.x_Mp2_OX;COM.Stk1.Tank2_OX.x_tank2_OX;COM.Stk1.Pipes2.x_pipes2;COM.Stk1.Aft_Struct2.x_aft_struct2;COM.Stk1.Engine2.x_engine2;COM.Stk1.Interstage.x_interstage;COM.Stk1.Prop1_FU.x_Mp1_FU;COM.Stk1.Tank1_FU.x_tank1_FU;COM.Stk1.Pipes1.x_pipes1;COM.Stk1.Prop1_OX.x_Mp1_OX;COM.Stk1.Tank1_OX.x_tank1_OX;COM.Stk1.Aft_Struct1.x_aft_struct1;COM.Stk1.Engine1.x_engine1;COM.Stk1.Fins.x_fins];
COM_vec_x0 = [COM.Stk1.Pay.x_pay;COM.Stk1.Fair.x_fair;COM.Stk1.Prac2.x_parachute_2;COM.Stk1.Prop2_FU.x_Mp2_FU0;COM.Stk1.Tank2_FU.x_tank2_FU;COM.Stk1.Prop2_OX.x_Mp2_OX0;COM.Stk1.Tank2_OX.x_tank2_OX;COM.Stk1.Pipes2.x_pipes2;COM.Stk1.Aft_Struct2.x_aft_struct2;COM.Stk1.Engine2.x_engine2;COM.Stk1.Interstage.x_interstage;COM.Stk1.Prop1_FU.x_Mp1_FU0;COM.Stk1.Tank1_FU.x_tank1_FU;COM.Stk1.Pipes1.x_pipes1;COM.Stk1.Prop1_OX.x_Mp1_OX0;COM.Stk1.Tank1_OX.x_tank1_OX;COM.Stk1.Aft_Struct1.x_aft_struct1;COM.Stk1.Engine1.x_engine1;COM.Stk1.Fins.x_fins];
%COM_vec_x0 = [COM.Pay.x_pay;COM.Fair.x_fair;COM.Prac2.x_parachute_2;COM.Prop2_FU.x_Mp2_FU0;COM.Tank2_FU.x_tank2_FU;COM.Prop2_OX.x_Mp2_OX0;COM.Tank2_OX.x_tank2_OX;COM.Pipes2.x_pipes2;COM.Aft_Struct2.x_aft_struct2;COM.Engine2.x_engine2;COM.Interstage.x_interstage;COM.Prop1_FU.x_Mp1_FU0;COM.Tank1_FU.x_tank1_FU;COM.Pipes1.x_pipes1;COM.Prop1_OX.x_Mp1_OX0;COM.Tank1_OX.x_tank1_OX;COM.Aft_Struct1.x_aft_struct1;COM.Engine1.x_engine1;COM.Fins.x_fins];

% No residuals of propellant in tanks

if COM.Stk1.Prop2_FU.m_prop2_FU <= 0

COM_vec_m(4,1)=0;
COM_vec_x(4,1)=0;
COM.Stk1.Prop2_FU.m_prop2_FU = 0;
COM.Stk1.Prop2_FU.x_Mp2_FU = 0;
end

if COM.Stk1.Prop2_OX.m_prop2_OX <= 0

COM_vec_m(6,1)=0;
COM_vec_x(6,1)=0;
COM.Stk1.Prop2_OX.m_prop2_OX = 0;
COM.Stk1.Prop2_OX.x_Mp2_OX = 0;
end

if COM.Stk1.Prop1_FU.m_prop1_FU <= 0

COM_vec_m(12,1)=0;
COM_vec_x(12,1)=0;
COM.Stk1.Prop1_FU.m_prop1_FU = 0;
COM.Stk1.Prop1_FU.x_Mp1_FU = 0;

end

if COM.Stk1.Prop1_OX.m_prop1_OX <= 0

COM_vec_m(15,1)=0;
COM_vec_x(15,1)=0;
COM.Stk1.Prop1_OX.m_prop1_OX = 0;
COM.Stk1.Prop1_OX.x_Mp1_OX = 0;


end

COM.Stk1.TR.M_LV = (sum(COM_vec_m));
COM.Stk1.TR.x_LV = (dot(COM_vec_m,COM_vec_x))/(sum(COM_vec_m));

% C.O.M. position of LV at t = t1 + t2 and at t = 0 (wrt nose)

COM.Stk1.TR.x_LV = (dot(COM_vec_m,COM_vec_x))/(sum(COM_vec_m)); % [m] t = t1 + t2
COM.Stk1.TR.x_LV0 = (dot(COM_vec_m0,COM_vec_x0))/(sum(COM_vec_m0)); % [m] t = 0

% Global mass of LV at t = t1 + t2 and at t = 0 (wrt nose)

COM.Stk1.TR.M_LV = (sum(COM_vec_m)); % [kg] t = t1 + t2
COM.Stk1.TR.M_LV0 = (sum(COM_vec_m0));% [kg] t = 0

% Propellant mass of LV at t = t1 + t2 and at t = 0 (wrt nose)

COM.Stk1.TR.M_prop = COM_vec_m(4,1) +COM_vec_m(6,1)+COM_vec_m(12,1) + COM_vec_m(15,1); % [kg] t = t1 + t2
COM.Stk1.TR.M_prop0 = COM.Stk1.Prop2_FU.m_prop2_FU0 + COM.Stk1.Prop2_OX.m_prop2_OX0 + COM.Stk1.Prop1_FU.m_prop1_FU0 +COM.Stk1.Prop1_OX.m_prop1_OX0;% [kg] t = 0

%% Moments of Inertia

% First the distance wrt to the C.O.M. is computed: (0: at t=0, [-]: at t)
% [m]

MoI.Stk1.Pay.x_pay = abs(COM.Stk1.TR.x_LV - COM.Stk1.Pay.x_pay);
MoI.Stk1.Pay.x_pay0 = abs(COM.Stk1.TR.x_LV0 - COM.Stk1.Pay.x_pay);

MoI.Stk1.Fair.x_fair = abs(COM.Stk1.TR.x_LV - COM.Stk1.Fair.x_fair);
MoI.Stk1.Fair.x_fair0 = abs(COM.Stk1.TR.x_LV0 - COM.Stk1.Fair.x_fair);

MoI.Stk1.Prac2.x_parachute_2 = abs(COM.Stk1.TR.x_LV - COM.Stk1.Prac2.x_parachute_2);
MoI.Stk1.Prac2.x_parachute_20 = abs(COM.Stk1.TR.x_LV0 - COM.Stk1.Prac2.x_parachute_2);

MoI.Stk1.Prop2_FU.x_Mp2_FU = abs(COM.Stk1.TR.x_LV - COM.Stk1.Prop2_FU.x_Mp2_FU);
MoI.Stk1.Prop2_FU.x_Mp2_FU0 = abs(COM.Stk1.TR.x_LV0 - COM.Stk1.Prop2_FU.x_Mp2_FU0);

MoI.Stk1.Tank2_FU.x_tank2_FU= abs(COM.Stk1.TR.x_LV - COM.Stk1.Tank2_FU.x_tank2_FU);
MoI.Stk1.Tank2_FU.x_tank2_FU0= abs(COM.Stk1.TR.x_LV0 - COM.Stk1.Tank2_FU.x_tank2_FU);

MoI.Stk1.Prop2_OX.x_Mp2_OX = abs(COM.Stk1.TR.x_LV - COM.Stk1.Prop2_OX.x_Mp2_OX);
MoI.Stk1.Prop2_OX.x_Mp2_OX0 = abs(COM.Stk1.TR.x_LV0 - COM.Stk1.Prop2_OX.x_Mp2_OX0);

MoI.Stk1.Tank2_OX.x_tank2_OX= abs(COM.Stk1.TR.x_LV - COM.Stk1.Tank2_OX.x_tank2_OX);
MoI.Stk1.Tank2_OX.x_tank2_OX0= abs(COM.Stk1.TR.x_LV0 - COM.Stk1.Tank2_OX.x_tank2_OX);

MoI.Stk1.Pipes2.x_pipes2 = abs(COM.Stk1.TR.x_LV - COM.Stk1.Pipes2.x_pipes2);
MoI.Stk1.Pipes2.x_pipes20 = abs(COM.Stk1.TR.x_LV0 - COM.Stk1.Pipes2.x_pipes2);

MoI.Stk1.Aft_Struct2.x_aft_struct2= abs(COM.Stk1.TR.x_LV - COM.Stk1.Aft_Struct2.x_aft_struct2);
MoI.Stk1.Aft_Struct2.x_aft_struct20= abs(COM.Stk1.TR.x_LV0 - COM.Stk1.Aft_Struct2.x_aft_struct2);

MoI.Stk1.Engine2.x_engine2= abs(COM.Stk1.TR.x_LV - COM.Stk1.Engine2.x_engine2);
MoI.Stk1.Engine2.x_engine20= abs(COM.Stk1.TR.x_LV0 - COM.Stk1.Engine2.x_engine2);

MoI.Stk1.Interstage.x_interstage = abs(COM.Stk1.TR.x_LV - COM.Stk1.Interstage.x_interstage);
MoI.Stk1.Interstage.x_interstage0 = abs(COM.Stk1.TR.x_LV0 - COM.Stk1.Interstage.x_interstage);

MoI.Stk1.Prop1_FU.x_Mp1_FU = abs(COM.Stk1.TR.x_LV - COM.Stk1.Prop1_FU.x_Mp1_FU);
MoI.Stk1.Prop1_FU.x_Mp1_FU0 = abs(COM.Stk1.TR.x_LV0 - COM.Stk1.Prop1_FU.x_Mp1_FU0);

MoI.Stk1.Tank1_FU.x_tank1_FU= abs(COM.Stk1.TR.x_LV - COM.Stk1.Tank1_FU.x_tank1_FU);
MoI.Stk1.Tank1_FU.x_tank1_FU0= abs(COM.Stk1.TR.x_LV0 - COM.Stk1.Tank1_FU.x_tank1_FU);

MoI.Stk1.Pipes1.x_pipes1 = abs(COM.Stk1.TR.x_LV - COM.Stk1.Pipes1.x_pipes1);
MoI.Stk1.Pipes1.x_pipes10 = abs(COM.Stk1.TR.x_LV0 - COM.Stk1.Pipes1.x_pipes1);

MoI.Stk1.Prop1_OX.x_Mp1_OX = abs(COM.Stk1.TR.x_LV - COM.Stk1.Prop1_OX.x_Mp1_OX);
MoI.Stk1.Prop1_OX.x_Mp1_OX0 = abs(COM.Stk1.TR.x_LV0 - COM.Stk1.Prop1_OX.x_Mp1_OX0);

MoI.Stk1.Tank1_OX.x_tank1_OX= abs(COM.Stk1.TR.x_LV - COM.Stk1.Tank1_OX.x_tank1_OX);
MoI.Stk1.Tank1_OX.x_tank1_OX0= abs(COM.Stk1.TR.x_LV0 - COM.Stk1.Tank1_OX.x_tank1_OX);

MoI.Stk1.Aft_Struct1.x_aft_struct1= abs(COM.Stk1.TR.x_LV - COM.Stk1.Aft_Struct1.x_aft_struct1);
MoI.Stk1.Aft_Struct1.x_aft_struct10= abs(COM.Stk1.TR.x_LV0 - COM.Stk1.Aft_Struct1.x_aft_struct1);

MoI.Stk1.Engine1.x_engine1= abs(COM.Stk1.TR.x_LV - COM.Stk1.Engine1.x_engine1);
MoI.Stk1.Engine1.x_engine10= abs(COM.Stk1.TR.x_LV0 - COM.Stk1.Engine1.x_engine1);

MoI.Stk1.Fins.x_fins= abs(COM.Stk1.TR.x_LV - COM.Stk1.Fins.x_fins);
MoI.Stk1.Fins.x_fins0= abs(COM.Stk1.TR.x_LV0 - COM.Stk1.Fins.x_fins);

if COM.Stk1.Prop2_FU.m_prop2_FU<=0

MoI.Stk1.Prop2_FU.x_Mp2_FU = 0;

end

if COM.Stk1.Prop2_OX.m_prop2_OX<=0

MoI.Stk1.Prop2_OX.x_Mp2_OX =0;

end

if COM.Stk1.Prop1_FU.m_prop1_FU<=0

MoI.Stk1.Prop1_FU.x_Mp1_FU = 0;

end

if COM.Stk1.Prop1_OX.m_prop1_OX<=0 

MoI.Stk1.Prop1_OX.x_Mp1_OX =0;

end

% Inertia moments of components (J0) are computed only for the tanks, 
% interstage, fairing while all other quantities are treated as point masses:

% All in [kg/m^2]:

MoI.Stk1.Fair.J0_fair = (3/10)*TR.Fair.M_fair*((Diam_2)^2);

if strcmp(Tank2.OX.FLAG, 'Cyl') && strcmp(Tank2.FU.FLAG, 'Cyl')

t_2O = (Tank2.OX.t_lox)/(Tank2.OX.R_cyl_lox + Tank2.OX.t_lox);
t_2F = (Tank2.FU.t_fuel)/(Tank2.FU.R_cyl_fuel + Tank2.FU.t_fuel);

MoI.Stk1.Tank2_OX.J0_tank2_OX = (Tank2.OX.M_tot_tank_lox*((Tank2.OX.R_cyl_lox)^2)*(1 - t_2O + ((t_2O^2)/2))) + ((4/5)*Tank2.OX.M_spherical_cap_lox*(Tank2.OX.R_cyl_lox));
MoI.Stk1.Tank2_FU.J0_tank2_FU = (Tank2.FU.M_tot_tank_fuel*((Tank2.FU.R_cyl_fuel)^2)*(1 - t_2F + ((t_2F^2)/2))) + ((4/5)*Tank2.FU.M_spherical_cap_fuel*(Tank2.FU.R_cyl_fuel));

elseif strcmp(Tank2.OX.FLAG, 'Sphere') && strcmp(Tank2.FU.FLAG, 'Sphere')

MoI.Stk1.Tank2_OX.J0_tank2_OX = (2/5)*(Tank2.OX.M_tot_tank_lox)*((((Tank2.OX.R_sphere_lox+ Tank2.OX.t_lox)^5)-((Tank2.OX.R_sphere_lox)^5))/(((Tank2.OX.R_sphere_lox+ Tank2.OX.t_lox)^3)-((Tank2.OX.R_sphere_lox)^3))); % [kg*m^2]

MoI.Stk1.Tank2_FU.J0_tank2_FU = (2/5)*(Tank2.FU.M_tot_tank_fuel)*((((Tank2.FU.R_sphere_fuel+ Tank2.FU.t_fuel)^5)-((Tank2.FU.R_sphere_fuel)^5))/(((Tank2.FU.R_sphere_fuel+ Tank2.FU.t_fuel)^3)-((Tank2.FU.R_sphere_fuel)^3)));% [kg*m^2]

end

t_E2 = Engine2.t_2/(Engine2.R_ext);
MoI.Stk1.Engine2.J0_engine2 = Engine2.M*((Engine2.R_ext)^2)*(1 - t_E2 + ((t_E2^2)/2));

Rad_1 = Diam_1/2;
Rad_2 = Diam_2/2;

if Rad_1 == Rad_2

MoI.Stk1.Interstage.J0_interstage = 0;

elseif Rad_1~=Rad_2

MoI.Stk1.Interstage.J0_interstage = COM.Stk1.Interstage.m_interstage*( (((Rad_1^2) + (Rad_2^2))/4) + (((L_interstage^2)/18)*(1 + ((2*Rad_1*Rad_2)/((Rad_1+Rad_2)^2))) ) );
end

if strcmp(Tank1.OX.FLAG, 'Cyl') && strcmp(Tank1.FU.FLAG, 'Cyl')

t_1O = (Tank1.OX.t_lox)/(Tank1.OX.R_cyl_lox + Tank1.OX.t_lox);
t_1F = (Tank1.FU.t_fuel)/(Tank1.FU.R_cyl_fuel + Tank1.FU.t_fuel);

MoI.Stk1.Tank1_OX.J0_tank1_OX = (Tank1.OX.M_tot_tank_lox*((Tank1.OX.R_cyl_lox)^2)*(1 - t_1O + ((t_1O^2)/2))) + ((4/5)*Tank1.OX.M_spherical_cap_lox*(Tank1.OX.R_cyl_lox));
MoI.Stk1.Tank1_FU.J0_tank1_FU = (Tank1.FU.M_tot_tank_fuel*((Tank1.FU.R_cyl_fuel)^2)*(1 - t_1F + ((t_1F^2)/2))) + ((4/5)*Tank1.FU.M_spherical_cap_fuel*(Tank1.FU.R_cyl_fuel));

elseif strcmp(Tank1.OX.FLAG, 'Sphere') && strcmp(Tank1.FU.FLAG, 'Sphere')

MoI.Stk1.Tank1_OX.J0_tank1_OX = (2/5)*(Tank1.OX.M_tot_tank_lox)*((((Tank1.OX.R_sphere_lox+Tank1.OX.t_lox)^5)-((Tank1.OX.R_sphere_lox)^5))/(((Tank1.OX.R_sphere_lox+Tank1.OX.t_lox)^3)-((Tank1.OX.R_sphere_lox)^3)));
MoI.Stk1.Tank1_FU.J0_tank1_FU = (2/5)*(Tank1.FU.M_tot_tank_fuel)*((((Tank1.FU.R_sphere_fuel+ Tank1.FU.t_fuel)^5)-((Tank1.FU.R_sphere_fuel)^5))/(((Tank1.FU.R_sphere_fuel+ Tank1.FU.t_fuel)^3)-((Tank1.FU.R_sphere_fuel)^3)));

end

t_E1 = Engine1.t_1/(Engine1.R_ext);
MoI.Stk1.Engine1.J0_engine1 = Engine1.M*((Engine1.R_ext)^2)*(1 - t_E1 + ((t_E1^2)/2));

% Vector assembly (careful with associating right components)

MoI_vec_J0 = [0;MoI.Stk1.Fair.J0_fair;0;0;MoI.Stk1.Tank2_FU.J0_tank2_FU;0;MoI.Stk1.Tank2_OX.J0_tank2_OX;0;0;MoI.Stk1.Engine2.J0_engine2;MoI.Stk1.Interstage.J0_interstage;0;MoI.Stk1.Tank1_FU.J0_tank1_FU;0;0;MoI.Stk1.Tank1_OX.J0_tank1_OX;0;MoI.Stk1.Engine1.J0_engine1;0];

MoI_vec_x = [MoI.Stk1.Pay.x_pay;MoI.Stk1.Fair.x_fair;MoI.Stk1.Prac2.x_parachute_2;MoI.Stk1.Prop2_FU.x_Mp2_FU;MoI.Stk1.Tank2_FU.x_tank2_FU;MoI.Stk1.Prop2_OX.x_Mp2_OX;MoI.Stk1.Tank2_OX.x_tank2_OX;MoI.Stk1.Pipes2.x_pipes2;MoI.Stk1.Aft_Struct2.x_aft_struct2;MoI.Stk1.Engine2.x_engine2;MoI.Stk1.Interstage.x_interstage;MoI.Stk1.Prop1_FU.x_Mp1_FU;MoI.Stk1.Tank1_FU.x_tank1_FU;MoI.Stk1.Pipes1.x_pipes1;MoI.Stk1.Prop1_OX.x_Mp1_OX;MoI.Stk1.Tank1_OX.x_tank1_OX;MoI.Stk1.Aft_Struct1.x_aft_struct1;MoI.Stk1.Engine1.x_engine1;MoI.Stk1.Fins.x_fins];
MoI_vec_x0 = [MoI.Stk1.Pay.x_pay0;MoI.Stk1.Fair.x_fair0;MoI.Stk1.Prac2.x_parachute_20;MoI.Stk1.Prop2_FU.x_Mp2_FU0;MoI.Stk1.Tank2_FU.x_tank2_FU0;MoI.Stk1.Prop2_OX.x_Mp2_OX0;MoI.Stk1.Tank2_OX.x_tank2_OX0;MoI.Stk1.Pipes2.x_pipes20;MoI.Stk1.Aft_Struct2.x_aft_struct20;MoI.Stk1.Engine2.x_engine20;MoI.Stk1.Interstage.x_interstage0;MoI.Stk1.Prop1_FU.x_Mp1_FU0;MoI.Stk1.Tank1_FU.x_tank1_FU0;MoI.Stk1.Pipes1.x_pipes10;MoI.Stk1.Prop1_OX.x_Mp1_OX0;MoI.Stk1.Tank1_OX.x_tank1_OX0;MoI.Stk1.Aft_Struct1.x_aft_struct10;MoI.Stk1.Engine1.x_engine10;MoI.Stk1.Fins.x_fins0];

MoI_vec_xquad = MoI_vec_x.*MoI_vec_x;
MoI_vec_xquad0 =MoI_vec_x0.*MoI_vec_x0;

MoI_vec_mxquad = dot(COM_vec_m,MoI_vec_xquad);
MoI_vec_mxquad0 = dot(COM_vec_m0,MoI_vec_xquad0);

% Moments of inertia of our LL at t = t1 + t2 and at t = 0

MoI.Stk1.TR.Jy = MoI_vec_mxquad + sum(MoI_vec_J0); % [kg m^2] t = t1 + t2

MoI.Stk1.TR.Jy0 = MoI_vec_mxquad0+ sum(MoI_vec_J0); % [kg m^2] t = 0

%% STACK 2:

% STACK 1 has been jettinsoned, so the COM and MoI are recomputed:
% x = 0: tip of nose

elseif t.t1>=t.t_burn1 

Ms2_stage_ratio = TR.M.Ms2/TR.M.Ms1;

M_avionics_s2 = M_avionics*Ms2_stage_ratio;
M_cables_s2 = M_cables*Ms2_stage_ratio;

t2 = t.t2;
t1 = t.t_burn1;

x_fair = TR.Fair.fair_length*(1+0.05);

if strcmp(Tank2.OX.FLAG, 'Cyl') && strcmp(Tank2.FU.FLAG, 'Cyl')
% Top FU, bottom OX
L_fwrd_skirt_s2 = ((1/3)*Diam_2) + Tank2.FU.H_dome;
L_cyl_tank_FU2 = Tank2.FU.H_cyl;
L_connection_2 = ((1/4)*Diam_2) + (Tank2.OX.H_dome + Tank2.FU.H_dome);
L_cyl_tank_OX2 = Tank2.OX.H_cyl;
L_dwrd_skirt_s2 = ((1/3)*Diam_2) + Tank2.OX.H_dome;

x_fwrd_skirt_s2 = x_fair + L_fwrd_skirt_s2;
x_cyl_22 = x_fwrd_skirt_s2 + L_cyl_tank_FU2;
x_connection_2 = x_cyl_22 + L_connection_2;
x_cyl_21 = x_connection_2 + L_cyl_tank_OX2;
x_dwrd_skirt_s2 = x_cyl_21 + L_dwrd_skirt_s2;

elseif strcmp(Tank2.OX.FLAG, 'Sphere') && strcmp(Tank2.FU.FLAG, 'Sphere')

   L_fwrd_skirt_s2 = (Tank2.FU.R_sphere_fuel/2)*(1+0.05);
   L_rem_sphere_FU2 = Tank2.FU.R_sphere_fuel;
   L_connection_2 = ((1/4)*Diam_2) + ((Tank2.OX.R_sphere_lox/2) + (Tank2.FU.R_sphere_fuel/2));
   L_rem_sphere_OX2 = Tank2.OX.R_sphere_lox;
   L_dwrd_skirt_s2 = (Tank2.OX.R_sphere_lox/2) + ((1/3)*Diam_2);

x_fwrd_skirt_s2 = x_fair + L_fwrd_skirt_s2;
x_sph_22 = x_fwrd_skirt_s2 + L_rem_sphere_FU2;
x_connection_2 = x_sph_22 + L_connection_2;
x_sph_21 = x_connection_2 + L_rem_sphere_OX2;
x_dwrd_skirt_s2 = x_sph_21 + L_dwrd_skirt_s2;

end

x_end_2 = x_dwrd_skirt_s2 + L_nozzle2;  % new end of LV, 2^nd stack end

%% C.O.M. positioning for each component

% x=0 at tip

COM.Stk2.Pay.x_pay = TR.Fair.fair_length/2;

COM.Stk2.Fair.x_fair = TR.Fair.fair_length*(2/3);

COM.Stk2.Prac2.x_parachute_2 = x_fair + L_fwrd_skirt_s2/2;

if strcmp(Tank2.OX.FLAG, 'Cyl') && strcmp(Tank2.FU.FLAG, 'Cyl')

COM.Stk2.Prop2_FU.x_Mp2_FU0 = x_fwrd_skirt_s2+ L_cyl_tank_FU2/2;
COM.Stk2.Prop2_FU.x_Mp2_FU = COM.Stk2.Prop2_FU.x_Mp2_FU0 - (((m_dotfuel_2*t2)/(Tank2.FU.M_fuel))*COM.Stk2.Prop2_FU.x_Mp2_FU0);
COM.Stk2.Tank2_FU.x_tank2_FU = COM.Stk2.Prop2_FU.x_Mp2_FU0;

COM.Stk2.Prop2_OX.x_Mp2_OX0 = x_connection_2 + L_cyl_tank_OX2/2;
COM.Stk2.Prop2_OX.x_Mp2_OX = COM.Stk2.Prop2_OX.x_Mp2_OX0 - (((m_dotlox_2*t2)/(Tank2.OX.M_lox))*COM.Stk2.Prop2_OX.x_Mp2_OX0);
COM.Stk2.Tank2_OX.x_tank2_OX = COM.Stk2.Prop2_OX.x_Mp2_OX0;

COM.Stk2.Pipes2.x_pipes2 = x_cyl_22 + Tank2.FU.H_dome; % about interface between 2 tanks

COM.Stk2.Aft_Struct2.x_aft_struct2 = x_cyl_21 + L_dwrd_skirt_s2/2;

COM.Stk2.Engine2.x_engine2 = x_cyl_21 + ((L_dwrd_skirt_s2)*(2/3));

elseif strcmp(Tank2.OX.FLAG, 'Sphere') && strcmp(Tank2.FU.FLAG, 'Sphere')

COM.Stk2.Prop2_FU.x_Mp2_FU0 =  x_fwrd_skirt_s2+ L_rem_sphere_FU2/2;
COM.Stk2.Prop2_FU.x_Mp2_FU = COM.Stk2.Prop2_FU.x_Mp2_FU0 - (((m_dotfuel_2*t2)/(Tank2.FU.M_fuel))*COM.Stk2.Prop2_FU.x_Mp2_FU0);
COM.Stk2.Tank2_FU.x_tank2_FU = COM.Stk2.Prop2_FU.x_Mp2_FU0;

COM.Stk2.Prop2_OX.x_Mp2_OX0 = x_connection_2 + L_rem_sphere_OX2/2;
COM.Stk2.Prop2_OX.x_Mp2_OX = COM.Stk2.Prop2_OX.x_Mp2_OX0 - (((m_dotlox_2*t2)/(Tank2.OX.M_lox))*COM.Stk2.Prop2_OX.x_Mp2_OX0);
COM.Stk2.Tank2_OX.x_tank2_OX = COM.Stk2.Prop2_OX.x_Mp2_OX0;

COM.Stk2.Pipes2.x_pipes2 = x_sph_22 + Tank2.FU.R_sphere_fuel/2; % about interface between 2 tanks

COM.Stk2.Aft_Struct2.x_aft_struct2 = x_sph_21 + L_dwrd_skirt_s2/2;

COM.Stk2.Engine2.x_engine2 = x_sph_21 + ((L_dwrd_skirt_s2)*(2/3));

end

%% Mass association to every component of 2nd stack (as before, no change in their mass)

COM.Stk2.Pay.m_pay = 250;

COM.Stk2.Fair.m_fair = TR.Fair.M_fair+ M_cables_distributed;

COM.Stk2.Prac2.m_parachute_2 = 0.1*TR.M.Ms2+ M_cables_distributed + M_avionics_s2/2;

COM.Stk2.Tank2_FU.m_tank2_FU = Tank2.FU.M_tot_tank_fuel+ M_cables_distributed;

COM.Stk2.Tank2_OX.m_tank2_OX = Tank2.OX.M_tot_tank_lox+ M_cables_distributed;

COM.Stk2.Prop2_FU.m_prop2_FU0 = Tank2.FU.M_fuel;
COM.Stk2.Prop2_FU.m_prop2_FU = COM.Stk2.Prop2_FU.m_prop2_FU0 - (m_dotfuel_2*t2);

COM.Stk2.Prop2_OX.m_prop2_OX0 = Tank2.OX.M_lox;
COM.Stk2.Prop2_OX.m_prop2_OX = COM.Stk2.Prop2_OX.m_prop2_OX0 - (m_dotlox_2*t2);

COM.Stk2.Pipes2.m_pipes2 = 10+ M_cables_distributed;

COM.Stk2.Aft_Struct2.m_aft_struct2 = 100+ M_cables_distributed+ M_avionics_s2/2;

COM.Stk2.Engine2.m_engine2 = 300+ M_cables_distributed;

%% Computation of position of COM wrt nose and M_tb1 (at t = tb1)

% Vector assembly ( careful with associating right components)

COM_vec_m = [COM.Stk2.Pay.m_pay;COM.Stk2.Fair.m_fair;COM.Stk2.Prac2.m_parachute_2;COM.Stk2.Prop2_FU.m_prop2_FU;COM.Stk2.Tank2_FU.m_tank2_FU;COM.Stk2.Prop2_OX.m_prop2_OX;COM.Stk2.Tank2_OX.m_tank2_OX;COM.Stk2.Pipes2.m_pipes2;COM.Stk2.Aft_Struct2.m_aft_struct2;COM.Stk2.Engine2.m_engine2];
COM_vec_m_tb1 = [COM.Stk2.Pay.m_pay;COM.Stk2.Fair.m_fair;COM.Stk2.Prac2.m_parachute_2;COM.Stk2.Prop2_FU.m_prop2_FU0;COM.Stk2.Tank2_FU.m_tank2_FU;COM.Stk2.Prop2_OX.m_prop2_OX0;COM.Stk2.Tank2_OX.m_tank2_OX;COM.Stk2.Pipes2.m_pipes2;COM.Stk2.Aft_Struct2.m_aft_struct2;COM.Stk2.Engine2.m_engine2];
COM_vec_x = [COM.Stk2.Pay.x_pay;COM.Stk2.Fair.x_fair;COM.Stk2.Prac2.x_parachute_2;COM.Stk2.Prop2_FU.x_Mp2_FU;COM.Stk2.Tank2_FU.x_tank2_FU;COM.Stk2.Prop2_OX.x_Mp2_OX;COM.Stk2.Tank2_OX.x_tank2_OX;COM.Stk2.Pipes2.x_pipes2;COM.Stk2.Aft_Struct2.x_aft_struct2;COM.Stk2.Engine2.x_engine2];
COM_vec_x_tb1 = [COM.Stk2.Pay.x_pay;COM.Stk2.Fair.x_fair;COM.Stk2.Prac2.x_parachute_2;COM.Stk2.Prop2_FU.x_Mp2_FU0;COM.Stk2.Tank2_FU.x_tank2_FU;COM.Stk2.Prop2_OX.x_Mp2_OX0;COM.Stk2.Tank2_OX.x_tank2_OX;COM.Stk2.Pipes2.x_pipes2;COM.Stk2.Aft_Struct2.x_aft_struct2;COM.Stk2.Engine2.x_engine2];

% No residuals of propellant in tanks

if COM.Stk2.Prop2_FU.m_prop2_FU<=0

COM_vec_m(4,1)=0;
COM_vec_x(4,1)=0;
COM.Stk2.Prop2_FU.m_prop2_FU = 0;
COM.Stk2.Prop2_FU.x_Mp2_FU = 0;
end

if COM.Stk2.Prop2_OX.m_prop2_OX<=0

COM_vec_m(6,1)=0;
COM_vec_x(6,1)=0;
COM.Stk2.Prop2_OX.m_prop2_OX = 0;
COM.Stk2.Prop2_OX.x_Mp2_OX = 0;
end

COM.Stk2.TR.M_LV = (sum(COM_vec_m));
COM.Stk2.TR.x_LV = (dot(COM_vec_m,COM_vec_x))/(sum(COM_vec_m));

% C.O.M. position of LV at t = tb1 + t2 and at t = tb1 (wrt nose)

COM.Stk2.TR.x_LV = (dot(COM_vec_m,COM_vec_x))/(sum(COM_vec_m)); % [m] t = tb1 + t2
COM.Stk2.TR.x_LV_tb1 = (dot(COM_vec_m_tb1,COM_vec_x_tb1))/(sum(COM_vec_m_tb1)); % [m] t = tb1

% Global mass of LV at t = tb1 + t2 and at t = tb1 (wrt nose)

COM.Stk2.TR.M_LV = (sum(COM_vec_m)); % [kg] t = tb1 + t2
COM.Stk2.TR.M_LV_tb1 = (sum(COM_vec_m_tb1));% [kg] t = tb1

% Propellant mass of LV at t = tb1 + t2 and at t = tb1 (wrt nose)

COM.Stk2.TR.M_prop = COM_vec_m(4,1) +COM_vec_m(6,1); % [kg] t = tb1 + t2
COM.Stk2.TR.M_prop_tb1 = COM.Stk2.Prop2_FU.m_prop2_FU0 + COM.Stk2.Prop2_OX.m_prop2_OX0;% [kg] t = tb1

%% Moments of Inertia

% First the distance wrt to the C.O.M. is computed: (tb1: at t=tb1, [-]: at t)
% [m]

MoI.Stk2.Pay.x_pay = abs(COM.Stk2.TR.x_LV - COM.Stk2.Pay.x_pay);
MoI.Stk2.Pay.x_pay_tb1 = abs(COM.Stk2.TR.x_LV_tb1 - COM.Stk2.Pay.x_pay);

MoI.Stk2.Fair.x_fair = abs(COM.Stk2.TR.x_LV - COM.Stk2.Fair.x_fair);
MoI.Stk2.Fair.x_fair_tb1 = abs(COM.Stk2.TR.x_LV_tb1 - COM.Stk2.Fair.x_fair);

MoI.Stk2.Prac2.x_parachute_2 = abs(COM.Stk2.TR.x_LV - COM.Stk2.Prac2.x_parachute_2);
MoI.Stk2.Prac2.x_parachute_2_tb1 = abs(COM.Stk2.TR.x_LV_tb1 - COM.Stk2.Prac2.x_parachute_2);

MoI.Stk2.Prop2_FU.x_Mp2_FU = abs(COM.Stk2.TR.x_LV - COM.Stk2.Prop2_FU.x_Mp2_FU);
MoI.Stk2.Prop2_FU.x_Mp2_FU_tb1 = abs(COM.Stk2.TR.x_LV_tb1 - COM.Stk2.Prop2_FU.x_Mp2_FU0);

MoI.Stk2.Tank2_FU.x_tank2_FU= abs(COM.Stk2.TR.x_LV - COM.Stk2.Tank2_FU.x_tank2_FU);
MoI.Stk2.Tank2_FU.x_tank2_FU_tb1= abs(COM.Stk2.TR.x_LV_tb1 - COM.Stk2.Tank2_FU.x_tank2_FU);

MoI.Stk2.Prop2_OX.x_Mp2_OX = abs(COM.Stk2.TR.x_LV - COM.Stk2.Prop2_OX.x_Mp2_OX);
MoI.Stk2.Prop2_OX.x_Mp2_OX_tb1 = abs(COM.Stk2.TR.x_LV_tb1 - COM.Stk2.Prop2_OX.x_Mp2_OX0);

MoI.Stk2.Tank2_OX.x_tank2_OX= abs(COM.Stk2.TR.x_LV - COM.Stk2.Tank2_OX.x_tank2_OX);
MoI.Stk2.Tank2_OX.x_tank2_OX_tb1= abs(COM.Stk2.TR.x_LV_tb1 - COM.Stk2.Tank2_OX.x_tank2_OX);

MoI.Stk2.Pipes2.x_pipes2 = abs(COM.Stk2.TR.x_LV - COM.Stk2.Pipes2.x_pipes2);
MoI.Stk2.Pipes2.x_pipes2_tb1 = abs(COM.Stk2.TR.x_LV_tb1 - COM.Stk2.Pipes2.x_pipes2);

MoI.Stk2.Aft_Struct2.x_aft_struct2= abs(COM.Stk2.TR.x_LV - COM.Stk2.Aft_Struct2.x_aft_struct2);
MoI.Stk2.Aft_Struct2.x_aft_struct2_tb1= abs(COM.Stk2.TR.x_LV_tb1 - COM.Stk2.Aft_Struct2.x_aft_struct2);

MoI.Stk2.Engine2.x_engine2= abs(COM.Stk2.TR.x_LV - COM.Stk2.Engine2.x_engine2);
MoI.Stk2.Engine2.x_engine2_tb1= abs(COM.Stk2.TR.x_LV_tb1 - COM.Stk2.Engine2.x_engine2);

if COM.Stk2.Prop2_FU.m_prop2_FU<=0

MoI.Stk2.Prop2_FU.x_Mp2_FU = 0;

end

if COM.Stk2.Prop2_OX.m_prop2_OX<=0

MoI.Stk2.Prop2_OX.x_Mp2_OX =0;

end

% Inertia moments of components (J0) are computed only for the tanks, 
% interstage, fairing while all other quantities are treated as point masses:

% All in [kg/m^2]:

MoI.Stk2.Fair.J0_fair = (3/10)*TR.Fair.M_fair*((Diam_2)^2);

if strcmp(Tank2.OX.FLAG, 'Cyl') && strcmp(Tank2.FU.FLAG, 'Cyl')

t_2O = (Tank2.OX.t_lox)/(Tank2.OX.R_cyl_lox + Tank2.OX.t_lox);
t_2F = (Tank2.FU.t_fuel)/(Tank2.FU.R_cyl_fuel + Tank2.FU.t_fuel);

MoI.Stk2.Tank2_OX.J0_tank2_OX = (Tank2.OX.M_tot_tank_lox*((Tank2.OX.R_cyl_lox)^2)*(1 - t_2O + ((t_2O^2)/2))) + ((4/5)*Tank2.OX.M_spherical_cap_lox*(Tank2.OX.R_cyl_lox));
MoI.Stk2.Tank2_FU.J0_tank2_FU = (Tank2.FU.M_tot_tank_fuel*((Tank2.FU.R_cyl_fuel)^2)*(1 - t_2F + ((t_2F^2)/2))) + ((4/5)*Tank2.FU.M_spherical_cap_fuel*(Tank2.FU.R_cyl_fuel));

elseif strcmp(Tank2.OX.FLAG, 'Sphere') && strcmp(Tank2.FU.FLAG, 'Sphere')

MoI.Stk2.Tank2_OX.J0_tank2_OX = (2/5)*(Tank2.OX.M_tot_tank_lox)*((((Tank2.OX.R_sphere_lox+ Tank2.OX.t_lox)^5)-((Tank2.OX.R_sphere_lox)^5))/(((Tank2.OX.R_sphere_lox+ Tank2.OX.t_lox)^3)-((Tank2.OX.R_sphere_lox)^3))); % [kg*m^2]

MoI.Stk2.Tank2_FU.J0_tank2_FU = (2/5)*(Tank2.FU.M_tot_tank_fuel)*((((Tank2.FU.R_sphere_fuel+ Tank2.FU.t_fuel)^5)-((Tank2.FU.R_sphere_fuel)^5))/(((Tank2.FU.R_sphere_fuel+ Tank2.FU.t_fuel)^3)-((Tank2.FU.R_sphere_fuel)^3)));% [kg*m^2]

end

t_E2 = Engine2.t_2/(Engine2.R_ext);
MoI.Stk2.Engine2.J0_engine2 = Engine2.M*((Engine2.R_ext)^2)*(1 - t_E2 + ((t_E2^2)/2));

% Vector assembly (careful with associating right components)

MoI_vec_J0 = [0;MoI.Stk2.Fair.J0_fair;0;0;MoI.Stk2.Tank2_FU.J0_tank2_FU;0;MoI.Stk2.Tank2_OX.J0_tank2_OX;0;0;MoI.Stk2.Engine2.J0_engine2];

MoI_vec_x = [MoI.Stk2.Pay.x_pay;MoI.Stk2.Fair.x_fair;MoI.Stk2.Prac2.x_parachute_2;MoI.Stk2.Prop2_FU.x_Mp2_FU;MoI.Stk2.Tank2_FU.x_tank2_FU;MoI.Stk2.Prop2_OX.x_Mp2_OX;MoI.Stk2.Tank2_OX.x_tank2_OX;MoI.Stk2.Pipes2.x_pipes2;MoI.Stk2.Aft_Struct2.x_aft_struct2;MoI.Stk2.Engine2.x_engine2];
MoI_vec_x_tb1 = [MoI.Stk2.Pay.x_pay_tb1;MoI.Stk2.Fair.x_fair_tb1;MoI.Stk2.Prac2.x_parachute_2_tb1;MoI.Stk2.Prop2_FU.x_Mp2_FU_tb1;MoI.Stk2.Tank2_FU.x_tank2_FU_tb1;MoI.Stk2.Prop2_OX.x_Mp2_OX_tb1;MoI.Stk2.Tank2_OX.x_tank2_OX_tb1;MoI.Stk2.Pipes2.x_pipes2_tb1;MoI.Stk2.Aft_Struct2.x_aft_struct2_tb1;MoI.Stk2.Engine2.x_engine2_tb1];

MoI_vec_xquad = MoI_vec_x.*MoI_vec_x;
MoI_vec_xquad_tb1 =MoI_vec_x_tb1.*MoI_vec_x_tb1;

MoI_vec_mxquad = dot(COM_vec_m,MoI_vec_xquad);
MoI_vec_mxquad_tb1 = dot(COM_vec_m_tb1,MoI_vec_xquad_tb1);

% Moments of inertia of our LL at t = tb1 + t2 and at t = tb1

MoI.Stk2.TR.Jy = MoI_vec_mxquad + sum(MoI_vec_J0); % [kg m^2] t = tb1 + t2

MoI.Stk2.TR.Jy0 = MoI_vec_mxquad_tb1+ sum(MoI_vec_J0); % [kg m^2] t = tb1

end 
end