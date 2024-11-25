function [X_COM,J_y] = COM_MoI_Stk1(M_prop_t,MASS)

OF = MASS.OF; %[-] Ox/Fu ratio for LOX-RP1

M_LOX = M_prop_t * OF / (1+OF);%[kg] mass of lox
M_FU = M_prop_t * 1  / (1+OF);%[kg] mass of rp1

%% 

L_nozzle2 = MASS.Engine2.L_nozzle;
L_nozzle1 = MASS.Engine1.L_nozzle;

Diam_1 = MASS.Diam_1;
Diam_2 = MASS.Diam_2;

x_fair = MASS.TR.Fair.fair_length;

if strcmp(MASS.Tank2.OX.FLAG, 'Cyl') && strcmp(MASS.Tank2.FU.FLAG, 'Cyl')
% Top FU, bottom OX
L_fwrd_skirt_s2 = (((1/3)*Diam_2) + MASS.Tank2.FU.H_dome);
L_cyl_tank_FU2 =  MASS.Tank2.FU.H_cyl;
L_connection_2 = ((1/4)*Diam_2) + (MASS.Tank2.OX.H_dome + MASS.Tank2.FU.H_dome);
%L_connection_2 =0;
L_cyl_tank_OX2 =  MASS.Tank2.OX.H_cyl;
L_dwrd_skirt_s2 = (((1/3)*Diam_2) +  MASS.Tank2.OX.H_dome);

x_fwrd_skirt_s2 = x_fair + L_fwrd_skirt_s2;
x_cyl_22 = x_fwrd_skirt_s2 + L_cyl_tank_FU2;
x_connection_2 = x_cyl_22 + L_connection_2;
x_cyl_21 = x_connection_2 + L_cyl_tank_OX2;
x_dwrd_skirt_s2 = x_cyl_21 + L_dwrd_skirt_s2;

elseif strcmp( MASS.Tank2.OX.FLAG, 'Sphere') && strcmp( MASS.Tank2.FU.FLAG, 'Sphere')

   L_fwrd_skirt_s2 = ( MASS.Tank2.FU.R_sphere_fuel/2);
   L_rem_sphere_FU2 =  MASS.Tank2.FU.R_sphere_fuel;
   %L_connection_2 =0;
   L_connection_2 = ((1/4)*Diam_2) + (( MASS.Tank2.OX.R_sphere_lox/2) + ( MASS.Tank2.FU.R_sphere_fuel/2));
   L_rem_sphere_OX2 =  MASS.Tank2.OX.R_sphere_lox;
   L_dwrd_skirt_s2 = (( MASS.Tank2.OX.R_sphere_lox/2) + ((1/3)*Diam_2));

x_fwrd_skirt_s2 = x_fair + L_fwrd_skirt_s2;
x_sph_22 = x_fwrd_skirt_s2 + L_rem_sphere_FU2;
x_connection_2 = x_sph_22 + L_connection_2;
x_sph_21 = x_connection_2 + L_rem_sphere_OX2;
x_dwrd_skirt_s2 = x_sph_21 + L_dwrd_skirt_s2;


end

L_interstage = Diam_1;%L_nozzle2*(1+0.1);

x_interstage = x_dwrd_skirt_s2 + L_interstage;

N_parts_s2 = 10; % 9 + half interstage
N_parts_s1 = 10; % 9 + half interstage

M_cables_s2_tot = MASS.M_cables_distributed*(x_dwrd_skirt_s2 + (L_interstage/2));
M_cables_s1_tot = MASS.M_cables_distributed*( MASS.TR.Length-(x_dwrd_skirt_s2 + (L_interstage/2)));

M_cables_s1 = M_cables_s1_tot/N_parts_s1;
M_cables_s2 = M_cables_s2_tot/N_parts_s2;


% 1 stage

if strcmp(MASS.Tank1.OX.FLAG, 'Cyl') && strcmp(MASS.Tank1.FU.FLAG, 'Cyl')
% Top FU, bottom OX
% Intarstage includes also dome of tank1 Fuel
L_cyl_tank_FU1 = MASS.Tank1.FU.H_cyl;
L_connection_1 = (((1/4)*Diam_1) + (MASS.Tank1.OX.H_dome + MASS.Tank1.FU.H_dome));
%L_connection_1 =0;
L_cyl_tank_OX1 = MASS.Tank1.OX.H_cyl;
L_end = Diam_1 + L_nozzle1;

x_cyl_12 = x_interstage + L_cyl_tank_FU1;
x_connection_1 = x_cyl_12 + L_connection_1;
x_cyl_11 = L_cyl_tank_OX1 + x_connection_1;
x_end = x_cyl_11 + L_end;

elseif strcmp(MASS.Tank1.OX.FLAG, 'Sphere') && strcmp(MASS.Tank1.FU.FLAG, 'Sphere')

   L_rem_sphere_FU1 = (Tank1.FU.R_sphere_fuel);
   %L_connection_1 =0;
   L_connection_1 = (((1/4)*Diam_1) + ((MASS.Tank1.OX.R_sphere_lox/2) + (MASS.Tank1.FU.R_sphere_fuel/2)));
   L_rem_sphere_OX1 = MASS.Tank1.OX.R_sphere_lox;
   L_end = L_nozzle1; 

x_sph_12 = x_interstage + L_rem_sphere_FU1;
x_connection_1 = x_sph_12 + L_connection_1;
x_sph_11 = L_rem_sphere_OX1 + x_connection_1;
x_end = x_sph_11 + L_end;

end

if x_end > MASS.TR.Length

x_end = MASS.TR.Length;

if strcmp(MASS.Tank1.OX.FLAG, 'Cyl') && strcmp(MASS.Tank1.FU.FLAG, 'Cyl') % correct end length

    L_end = MASS.TR.Length - x_cyl_11;

elseif strcmp(MASS.Tank1.OX.FLAG, 'Sphere') && strcmp(MASS.Tank1.FU.FLAG, 'Sphere')
    
L_end = MASS.TR.Length - x_sph_11;

end

end


%% C.O.M. positioning for each component

% x=0 at tip

COM.Stk1.Pay.x_pay = MASS.TR.Fair.fair_length/2;

COM.Stk1.Fair.x_fair = MASS.TR.Fair.fair_length*(2/3);

COM.Stk1.Frw_struct2.x_parachute_2 = x_fair + L_fwrd_skirt_s2/2;

if strcmp(MASS.Tank2.OX.FLAG, 'Cyl') && strcmp(MASS.Tank2.FU.FLAG, 'Cyl')

COM.Stk1.Prop2_FU.x_Mp2_FU0 = x_fwrd_skirt_s2+ L_cyl_tank_FU2/2;
COM.Stk1.Prop2_FU.x_Mp2_FU = COM.Stk1.Prop2_FU.x_Mp2_FU0;
COM.Stk1.Tank2_FU.x_tank2_FU = COM.Stk1.Prop2_FU.x_Mp2_FU0;

COM.Stk1.Prop2_OX.x_Mp2_OX0 = x_connection_2 + L_cyl_tank_OX2/2;
COM.Stk1.Prop2_OX.x_Mp2_OX = COM.Stk1.Prop2_OX.x_Mp2_OX0;
COM.Stk1.Tank2_OX.x_tank2_OX = COM.Stk1.Prop2_OX.x_Mp2_OX0;

COM.Stk1.Pipes2.x_pipes2 = x_cyl_22 + MASS.Tank2.FU.H_dome; % about interface between 2 tanks

COM.Stk1.Aft_Struct2.x_aft_struct2 = x_cyl_21 + L_dwrd_skirt_s2/2;

COM.Stk1.Engine2.x_engine2 = x_cyl_21 + ((L_dwrd_skirt_s2)*(2/3));

M_insulation_tank2_FU = MASS.M_insulation_tank2_FU;
M_insulation_tank2_OX = MASS.M_insulation_tank2_OX;

elseif strcmp(MASS.Tank2.OX.FLAG, 'Sphere') && strcmp(MASS.Tank2.FU.FLAG, 'Sphere')

COM.Stk1.Prop2_FU.x_Mp2_FU0 =  x_fwrd_skirt_s2+ L_rem_sphere_FU2/2;
COM.Stk1.Prop2_FU.x_Mp2_FU = COM.Stk1.Prop2_FU.x_Mp2_FU0;
COM.Stk1.Tank2_FU.x_tank2_FU = COM.Stk1.Prop2_FU.x_Mp2_FU0;

COM.Stk1.Prop2_OX.x_Mp2_OX0 = x_connection_2 + L_rem_sphere_OX2/2;
COM.Stk1.Prop2_OX.x_Mp2_OX = COM.Stk1.Prop2_OX.x_Mp2_OX0;
COM.Stk1.Tank2_OX.x_tank2_OX = COM.Stk1.Prop2_OX.x_Mp2_OX0;

COM.Stk1.Pipes2.x_pipes2 = x_sph_22 + MASS.Tank2.FU.R_sphere_fuel/2; % about interface between 2 tanks

COM.Stk1.Aft_Struct2.x_aft_struct2 = x_sph_21 + L_dwrd_skirt_s2/2;

COM.Stk1.Engine2.x_engine2 = x_sph_21 + ((L_dwrd_skirt_s2)*(2/3));

M_insulation_tank2_FU = MASS.M_insulation_tank2_FU;
M_insulation_tank2_OX = MASS.M_insulation_tank2_OX;

end

if Diam_1==Diam_2

COM.Stk1.Interstage.x_interstage = x_dwrd_skirt_s2 + (L_interstage*(1/2)); % Position that accounts for mass of nozzle + parachute + part of avionics; 1/2: to be changed maybe

elseif Diam_1~=Diam_2
% Parameters of the conical trunk
R1 = Diam_1/2;    % Radius of the bottom base
R2 = Diam_2/2;    % Radius of the top base
h = L_interstage;    % Height of the frustum

% Calculate the center of mass along the vertical axis (z)
z_cm = (h / 4) * ((R1^2 + 2*R1*R2 + 3*R2^2) / (R1^2 + R1*R2 + R2^2));


COM.Stk1.Interstage.x_interstage = x_dwrd_skirt_s2 + z_cm;    

end

if strcmp(MASS.Tank1.OX.FLAG, 'Cyl') && strcmp(MASS.Tank1.FU.FLAG, 'Cyl')

COM.Stk1.Prop1_FU.x_Mp1_FU0 = x_interstage + MASS.Tank1.FU.H_cyl/2;
COM.Stk1.Prop1_FU.x_Mp1_FU = COM.Stk1.Prop1_FU.x_Mp1_FU0 - (((abs(M_FU - MASS.Tank1.FU.M_fuel))/(MASS.Tank1.FU.M_fuel))*COM.Stk1.Prop1_FU.x_Mp1_FU0);
COM.Stk1.Tank1_FU.x_tank1_FU = COM.Stk1.Prop1_FU.x_Mp1_FU0;

COM.Stk1.Pipes1.x_pipes1 = x_cyl_12 + MASS.Tank1.FU.H_dome;

COM.Stk1.Prop1_OX.x_Mp1_OX0 = x_connection_1 + MASS.Tank1.OX.H_cyl/2;
COM.Stk1.Prop1_OX.x_Mp1_OX = COM.Stk1.Prop1_OX.x_Mp1_OX0 - (((abs(M_LOX-MASS.Tank1.OX.M_lox))/(MASS.Tank1.OX.M_lox))*COM.Stk1.Prop1_OX.x_Mp1_OX0);
COM.Stk1.Tank1_OX.x_tank1_OX=COM.Stk1.Prop1_OX.x_Mp1_OX0;

M_insulation_tank1_FU = MASS.M_insulation_tank1_FU;
M_insulation_tank1_OX =  MASS.M_insulation_tank1_OX;


elseif strcmp( MASS.Tank1.OX.FLAG, 'Sphere') && strcmp( MASS.Tank1.FU.FLAG, 'Sphere')

COM.Stk1.Prop1_FU.x_Mp1_FU0 = x_interstage + L_rem_sphere_FU1/2;
COM.Stk1.Prop1_FU.x_Mp1_FU = M_FU;
COM.Stk1.Tank1_FU.x_tank1_FU = COM.Stk1.Prop1_FU.x_Mp1_FU0;


COM.Stk1.Pipes1.x_pipes1 = x_sph_12 + Tank1.FU.R_sphere_fuel/2;

COM.Stk1.Prop1_OX.x_Mp1_OX0 = x_connection_1 + L_rem_sphere_OX1/2;
COM.Stk1.Prop1_OX.x_Mp1_OX = M_LOX;
COM.Stk1.Tank1_OX.x_tank1_OX=COM.Stk1.Prop1_OX.x_Mp1_OX0;

M_insulation_tank1_FU = 2.88* (4*pi*Tank1.FU.R_sphere_fuel^2);
M_insulation_tank1_OX = 1.123*(4*pi*Tank1.OX.R_sphere_lox^2);


end

COM.Stk1.Aft_Struct1.x_aft_struct1 = x_cyl_11 + (L_end*(1/4)); % Sustaining structures

COM.Stk1.Engine1.x_engine1 = x_cyl_11 + (L_end*(1/2)); % Engine + nozzle

COM.Stk1.Fins.x_fins = x_cyl_11 + (L_end*(3/4)); % x coord of fins

%% Mass association to every quantity in COM struct:

COM.Stk1.Pay.m_pay = 250;

COM.Stk1.Fair.m_fair = MASS.TR.Fair.M_fair +M_cables_s2;

COM.Stk1.Frw_struct2.m_parachute_2 = 0.07*MASS.TR.M.Ms2+ M_cables_s2 + MASS.M_avionics_s2/2 + MASS.M_struct2/2;

COM.Stk1.Tank2_FU.m_tank2_FU = MASS.Tank2.FU.M_tot_tank_fuel+ M_cables_s2 + M_insulation_tank2_FU;

COM.Stk1.Tank2_OX.m_tank2_OX = MASS.Tank2.OX.M_tot_tank_lox+ M_cables_s2 + M_insulation_tank2_OX;

COM.Stk1.Prop2_FU.m_prop2_FU0 = MASS.Tank2.FU.M_fuel;
COM.Stk1.Prop2_FU.m_prop2_FU = COM.Stk1.Prop2_FU.m_prop2_FU0;

COM.Stk1.Prop2_OX.m_prop2_OX0 = MASS.Tank2.OX.M_lox;
COM.Stk1.Prop2_OX.m_prop2_OX = COM.Stk1.Prop2_OX.m_prop2_OX0;

COM.Stk1.Pipes2.m_pipes2 = 10+ M_cables_s2;

COM.Stk1.Aft_Struct2.m_aft_struct2 = MASS.M_struct2/2+ M_cables_s2+ MASS.M_avionics_s2/2;

COM.Stk1.Engine2.m_engine2 = MASS.Engine2.M+ M_cables_s2;

COM.Stk1.Interstage.m_interstage = 0.07*MASS.TR.M.Ms1 + (MASS.M_struct1+MASS.M_struct2)/2+ M_cables_s1+M_cables_s2+ MASS.M_avionics_s1/2;

COM.Stk1.Tank1_FU.m_tank1_FU = MASS.Tank1.FU.M_tot_tank_fuel+ M_cables_s1 + M_insulation_tank1_FU;

COM.Stk1.Tank1_OX.m_tank1_OX = MASS.Tank1.OX.M_tot_tank_lox+ M_cables_s1 + M_insulation_tank1_OX;

COM.Stk1.Prop1_FU.m_prop1_FU0 = MASS.Tank1.FU.M_fuel;
COM.Stk1.Prop1_FU.m_prop1_FU = M_FU;

COM.Stk1.Prop1_OX.m_prop1_OX0 = MASS.Tank1.OX.M_lox;
COM.Stk1.Prop1_OX.m_prop1_OX = M_LOX;
COM.Stk1.Pipes1.m_pipes1 = 10 + M_cables_s1;

COM.Stk1.Aft_Struct1.m_aft_struct1 = MASS.M_struct1+ M_cables_s1 +MASS.M_avionics_s1;

COM.Stk1.Engine1.m_engine1 = MASS.Engine1.M + M_cables_s1; % Engine + nozzle

COM.Stk1.Fins.m_fins = 10 + M_cables_s1; % x coord of fins

%% Computation of position of COM wrt nose and M0

% Vector assembly ( careful with associating right components)

COM_vec_m = [COM.Stk1.Pay.m_pay;COM.Stk1.Fair.m_fair;COM.Stk1.Frw_struct2.m_parachute_2;COM.Stk1.Prop2_FU.m_prop2_FU;COM.Stk1.Tank2_FU.m_tank2_FU;COM.Stk1.Prop2_OX.m_prop2_OX;COM.Stk1.Tank2_OX.m_tank2_OX;COM.Stk1.Pipes2.m_pipes2;COM.Stk1.Aft_Struct2.m_aft_struct2;COM.Stk1.Engine2.m_engine2;COM.Stk1.Interstage.m_interstage;COM.Stk1.Prop1_FU.m_prop1_FU;COM.Stk1.Tank1_FU.m_tank1_FU;COM.Stk1.Pipes1.m_pipes1;COM.Stk1.Prop1_OX.m_prop1_OX;COM.Stk1.Tank1_OX.m_tank1_OX;COM.Stk1.Aft_Struct1.m_aft_struct1;COM.Stk1.Engine1.m_engine1;COM.Stk1.Fins.m_fins];
COM_vec_m0 = [COM.Stk1.Pay.m_pay;COM.Stk1.Fair.m_fair;COM.Stk1.Frw_struct2.m_parachute_2;COM.Stk1.Prop2_FU.m_prop2_FU0;COM.Stk1.Tank2_FU.m_tank2_FU;COM.Stk1.Prop2_OX.m_prop2_OX0;COM.Stk1.Tank2_OX.m_tank2_OX;COM.Stk1.Pipes2.m_pipes2;COM.Stk1.Aft_Struct2.m_aft_struct2;COM.Stk1.Engine2.m_engine2;COM.Stk1.Interstage.m_interstage;COM.Stk1.Prop1_FU.m_prop1_FU0;COM.Stk1.Tank1_FU.m_tank1_FU;COM.Stk1.Pipes1.m_pipes1;COM.Stk1.Prop1_OX.m_prop1_OX0;COM.Stk1.Tank1_OX.m_tank1_OX;COM.Stk1.Aft_Struct1.m_aft_struct1;COM.Stk1.Engine1.m_engine1;COM.Stk1.Fins.m_fins];
%COM_vec_m0 = [COM.Pay.m_pay;COM.Fair.m_fair;COM.Frw_struct2.m_parachute_2;COM.Prop2_FU.m_prop2_FU0;COM.Tank2_FU.m_tank2_FU;COM.Prop2_OX.m_prop2_OX0;COM.Tank2_OX.m_tank2_OX;COM.Pipes2.m_pipes2;COM.Aft_Struct2.m_aft_struct2;COM.Engine2.m_engine2;COM.Interstage.m_interstage;COM.Prop1_FU.m_prop1_FU0;COM.Tank1_FU.m_tank1_FU;COM.Pipes1.m_pipes1;COM.Prop1_OX.m_prop1_OX0;COM.Tank1_OX.m_tank1_OX;COM.Stk1.Aft_Struct1.m_aft_struct1;COM.Stk1.Engine1.m_engine1;COM.Stk1.Fins.m_fins];
COM_vec_x = [COM.Stk1.Pay.x_pay;COM.Stk1.Fair.x_fair;COM.Stk1.Frw_struct2.x_parachute_2;COM.Stk1.Prop2_FU.x_Mp2_FU;COM.Stk1.Tank2_FU.x_tank2_FU;COM.Stk1.Prop2_OX.x_Mp2_OX;COM.Stk1.Tank2_OX.x_tank2_OX;COM.Stk1.Pipes2.x_pipes2;COM.Stk1.Aft_Struct2.x_aft_struct2;COM.Stk1.Engine2.x_engine2;COM.Stk1.Interstage.x_interstage;COM.Stk1.Prop1_FU.x_Mp1_FU;COM.Stk1.Tank1_FU.x_tank1_FU;COM.Stk1.Pipes1.x_pipes1;COM.Stk1.Prop1_OX.x_Mp1_OX;COM.Stk1.Tank1_OX.x_tank1_OX;COM.Stk1.Aft_Struct1.x_aft_struct1;COM.Stk1.Engine1.x_engine1;COM.Stk1.Fins.x_fins];
COM_vec_x0 = [COM.Stk1.Pay.x_pay;COM.Stk1.Fair.x_fair;COM.Stk1.Frw_struct2.x_parachute_2;COM.Stk1.Prop2_FU.x_Mp2_FU0;COM.Stk1.Tank2_FU.x_tank2_FU;COM.Stk1.Prop2_OX.x_Mp2_OX0;COM.Stk1.Tank2_OX.x_tank2_OX;COM.Stk1.Pipes2.x_pipes2;COM.Stk1.Aft_Struct2.x_aft_struct2;COM.Stk1.Engine2.x_engine2;COM.Stk1.Interstage.x_interstage;COM.Stk1.Prop1_FU.x_Mp1_FU0;COM.Stk1.Tank1_FU.x_tank1_FU;COM.Stk1.Pipes1.x_pipes1;COM.Stk1.Prop1_OX.x_Mp1_OX0;COM.Stk1.Tank1_OX.x_tank1_OX;COM.Stk1.Aft_Struct1.x_aft_struct1;COM.Stk1.Engine1.x_engine1;COM.Stk1.Fins.x_fins];
%COM_vec_x0 = [COM.Pay.x_pay;COM.Fair.x_fair;COM.Frw_struct2.x_parachute_2;COM.Prop2_FU.x_Mp2_FU0;COM.Tank2_FU.x_tank2_FU;COM.Prop2_OX.x_Mp2_OX0;COM.Tank2_OX.x_tank2_OX;COM.Pipes2.x_pipes2;COM.Aft_Struct2.x_aft_struct2;COM.Engine2.x_engine2;COM.Interstage.x_interstage;COM.Prop1_FU.x_Mp1_FU0;COM.Tank1_FU.x_tank1_FU;COM.Pipes1.x_pipes1;COM.Prop1_OX.x_Mp1_OX0;COM.Tank1_OX.x_tank1_OX;COM.Aft_Struct1.x_aft_struct1;COM.Engine1.x_engine1;COM.Fins.x_fins];

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


X_COM = COM.Stk1.TR.x_LV;

%% Moments of Inertia

% First the distance wrt to the C.O.M. is computed: (0: at t=0, [-]: at t)
% [m]

MoI.Stk1.Pay.x_pay = abs(COM.Stk1.TR.x_LV - COM.Stk1.Pay.x_pay);

MoI.Stk1.Fair.x_fair = abs(COM.Stk1.TR.x_LV - COM.Stk1.Fair.x_fair);

MoI.Stk1.Frw_struct2.x_parachute_2 = abs(COM.Stk1.TR.x_LV - COM.Stk1.Frw_struct2.x_parachute_2);

MoI.Stk1.Prop2_FU.x_Mp2_FU = abs(COM.Stk1.TR.x_LV - COM.Stk1.Prop2_FU.x_Mp2_FU);

MoI.Stk1.Tank2_FU.x_tank2_FU= abs(COM.Stk1.TR.x_LV - COM.Stk1.Tank2_FU.x_tank2_FU);

MoI.Stk1.Prop2_OX.x_Mp2_OX = abs(COM.Stk1.TR.x_LV - COM.Stk1.Prop2_OX.x_Mp2_OX);

MoI.Stk1.Tank2_OX.x_tank2_OX= abs(COM.Stk1.TR.x_LV - COM.Stk1.Tank2_OX.x_tank2_OX);

MoI.Stk1.Pipes2.x_pipes2 = abs(COM.Stk1.TR.x_LV - COM.Stk1.Pipes2.x_pipes2);

MoI.Stk1.Aft_Struct2.x_aft_struct2= abs(COM.Stk1.TR.x_LV - COM.Stk1.Aft_Struct2.x_aft_struct2);

MoI.Stk1.Engine2.x_engine2= abs(COM.Stk1.TR.x_LV - COM.Stk1.Engine2.x_engine2);

MoI.Stk1.Interstage.x_interstage = abs(COM.Stk1.TR.x_LV - COM.Stk1.Interstage.x_interstage);

MoI.Stk1.Prop1_FU.x_Mp1_FU = abs(COM.Stk1.TR.x_LV - COM.Stk1.Prop1_FU.x_Mp1_FU);

MoI.Stk1.Tank1_FU.x_tank1_FU= abs(COM.Stk1.TR.x_LV - COM.Stk1.Tank1_FU.x_tank1_FU);

MoI.Stk1.Pipes1.x_pipes1 = abs(COM.Stk1.TR.x_LV - COM.Stk1.Pipes1.x_pipes1);

MoI.Stk1.Prop1_OX.x_Mp1_OX = abs(COM.Stk1.TR.x_LV - COM.Stk1.Prop1_OX.x_Mp1_OX);

MoI.Stk1.Tank1_OX.x_tank1_OX= abs(COM.Stk1.TR.x_LV - COM.Stk1.Tank1_OX.x_tank1_OX);

MoI.Stk1.Aft_Struct1.x_aft_struct1= abs(COM.Stk1.TR.x_LV - COM.Stk1.Aft_Struct1.x_aft_struct1);

MoI.Stk1.Engine1.x_engine1= abs(COM.Stk1.TR.x_LV - COM.Stk1.Engine1.x_engine1);

MoI.Stk1.Fins.x_fins= abs(COM.Stk1.TR.x_LV - COM.Stk1.Fins.x_fins);

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

MoI.Stk1.Fair.J0_fair = (3/10)*MASS.TR.Fair.M_fair*((Diam_2)^2);

if strcmp(MASS.Tank2.OX.FLAG, 'Cyl') && strcmp(MASS.Tank2.FU.FLAG, 'Cyl')

t_2O = (MASS.Tank2.OX.t_lox)/(MASS.Tank2.OX.R_cyl_lox + MASS.Tank2.OX.t_lox);
t_2F = (MASS.Tank2.FU.t_fuel)/(MASS.Tank2.FU.R_cyl_fuel + MASS.Tank2.FU.t_fuel);

MoI.Stk1.Tank2_OX.J0_tank2_OX = (MASS.Tank2.OX.M_tot_tank_lox*((MASS.Tank2.OX.R_cyl_lox)^2)*(1 - t_2O + ((t_2O^2)/2))) + ((4/5)*MASS.Tank2.OX.M_spherical_cap_lox*(MASS.Tank2.OX.R_cyl_lox));
MoI.Stk1.Tank2_FU.J0_tank2_FU = (MASS.Tank2.FU.M_tot_tank_fuel*((MASS.Tank2.FU.R_cyl_fuel)^2)*(1 - t_2F + ((t_2F^2)/2))) + ((4/5)*MASS.Tank2.FU.M_spherical_cap_fuel*(MASS.Tank2.FU.R_cyl_fuel));

elseif strcmp(MASS.Tank2.OX.FLAG, 'Sphere') && strcmp(MASS.Tank2.FU.FLAG, 'Sphere')

MoI.Stk1.Tank2_OX.J0_tank2_OX = (2/5)*(MASS.Tank2.OX.M_tot_tank_lox)*((((MASS.Tank2.OX.R_sphere_lox+ MASS.Tank2.OX.t_lox)^5)-((MASS.Tank2.OX.R_sphere_lox)^5))/(((MASS.Tank2.OX.R_sphere_lox+ MASS.Tank2.OX.t_lox)^3)-((MASS.Tank2.OX.R_sphere_lox)^3))); % [kg*m^2]

MoI.Stk1.Tank2_FU.J0_tank2_FU = (2/5)*(MASS.Tank2.FU.M_tot_tank_fuel)*((((MASS.Tank2.FU.R_sphere_fuel+ MASS.Tank2.FU.t_fuel)^5)-((MASS.Tank2.FU.R_sphere_fuel)^5))/(((MASS.Tank2.FU.R_sphere_fuel+ MASS.Tank2.FU.t_fuel)^3)-((MASS.Tank2.FU.R_sphere_fuel)^3)));% [kg*m^2]

end

Rad_1 = Diam_1/2;
Rad_2 = Diam_2/2;

if Rad_1 == Rad_2

MoI.Stk1.Interstage.J0_interstage = 0;

elseif Rad_1~=Rad_2

MoI.Stk1.Interstage.J0_interstage = COM.Stk1.Interstage.m_interstage*( (((Rad_1^2) + (Rad_2^2))/4) + (((L_interstage^2)/18)*(1 + ((2*Rad_1*Rad_2)/((Rad_1+Rad_2)^2))) ) );
end

if strcmp(MASS.Tank1.OX.FLAG, 'Cyl') && strcmp(MASS.Tank1.FU.FLAG, 'Cyl')

t_1O = (MASS.Tank1.OX.t_lox)/(MASS.Tank1.OX.R_cyl_lox + MASS.Tank1.OX.t_lox);
t_1F = (MASS.Tank1.FU.t_fuel)/(MASS.Tank1.FU.R_cyl_fuel + MASS.Tank1.FU.t_fuel);

MoI.Stk1.Tank1_OX.J0_tank1_OX = (MASS.Tank1.OX.M_tot_tank_lox*((MASS.Tank1.OX.R_cyl_lox)^2)*(1 - t_1O + ((t_1O^2)/2))) + ((4/5)*MASS.Tank1.OX.M_spherical_cap_lox*(MASS.Tank1.OX.R_cyl_lox));
MoI.Stk1.Tank1_FU.J0_tank1_FU = (MASS.Tank1.FU.M_tot_tank_fuel*((MASS.Tank1.FU.R_cyl_fuel)^2)*(1 - t_1F + ((t_1F^2)/2))) + ((4/5)*MASS.Tank1.FU.M_spherical_cap_fuel*(MASS.Tank1.FU.R_cyl_fuel));

elseif strcmp(MASS.Tank1.OX.FLAG, 'Sphere') && strcmp(MASS.Tank1.FU.FLAG, 'Sphere')

MoI.Stk1.Tank1_OX.J0_tank1_OX = (2/5)*(MASS.Tank1.OX.M_tot_tank_lox)*((((MASS.Tank1.OX.R_sphere_lox+MASS.Tank1.OX.t_lox)^5)-((MASS.Tank1.OX.R_sphere_lox)^5))/(((MASS.Tank1.OX.R_sphere_lox+MASS.Tank1.OX.t_lox)^3)-((MASS.Tank1.OX.R_sphere_lox)^3)));
MoI.Stk1.Tank1_FU.J0_tank1_FU = (2/5)*(MASS.Tank1.FU.M_tot_tank_fuel)*((((MASS.Tank1.FU.R_sphere_fuel+ MASS.Tank1.FU.t_fuel)^5)-((MASS.Tank1.FU.R_sphere_fuel)^5))/(((MASS.Tank1.FU.R_sphere_fuel+ MASS.Tank1.FU.t_fuel)^3)-((MASS.Tank1.FU.R_sphere_fuel)^3)));

end

% Vector assembly (careful with associating right components)

MoI_vec_J0 = [0;MoI.Stk1.Fair.J0_fair;0;0;MoI.Stk1.Tank2_FU.J0_tank2_FU;0;MoI.Stk1.Tank2_OX.J0_tank2_OX;0;0;0;MoI.Stk1.Interstage.J0_interstage;0;MoI.Stk1.Tank1_FU.J0_tank1_FU;0;0;MoI.Stk1.Tank1_OX.J0_tank1_OX;0;0;0];

MoI_vec_x = [MoI.Stk1.Pay.x_pay;MoI.Stk1.Fair.x_fair;MoI.Stk1.Frw_struct2.x_parachute_2;MoI.Stk1.Prop2_FU.x_Mp2_FU;MoI.Stk1.Tank2_FU.x_tank2_FU;MoI.Stk1.Prop2_OX.x_Mp2_OX;MoI.Stk1.Tank2_OX.x_tank2_OX;MoI.Stk1.Pipes2.x_pipes2;MoI.Stk1.Aft_Struct2.x_aft_struct2;MoI.Stk1.Engine2.x_engine2;MoI.Stk1.Interstage.x_interstage;MoI.Stk1.Prop1_FU.x_Mp1_FU;MoI.Stk1.Tank1_FU.x_tank1_FU;MoI.Stk1.Pipes1.x_pipes1;MoI.Stk1.Prop1_OX.x_Mp1_OX;MoI.Stk1.Tank1_OX.x_tank1_OX;MoI.Stk1.Aft_Struct1.x_aft_struct1;MoI.Stk1.Engine1.x_engine1;MoI.Stk1.Fins.x_fins];

MoI_vec_xquad = MoI_vec_x.*MoI_vec_x;

MoI_vec_mxquad = dot(COM_vec_m,MoI_vec_xquad);

% Moments of inertia of our LV at t = t1 + t2 and at t = 0

MoI.Stk1.TR.Jy = MoI_vec_mxquad + sum(MoI_vec_J0); % [kg m^2] t = t1 + t2

J_y = MoI.Stk1.TR.Jy;

end