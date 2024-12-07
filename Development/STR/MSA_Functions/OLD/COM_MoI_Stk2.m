function [X_COM,J_yaw,J_pitch,J_roll] = COM_MoI_Stk2(M_prop_t,MASS)

OF = MASS.OF; %[-] Ox/Fu ratio for LOX-RP1

M_LOX = M_prop_t * OF / (1+OF);%[kg] mass of lox
M_FU = M_prop_t * 1  / (1+OF);%[kg] mass of rp1

L_nozzle2 = MASS.Engine2.L_nozzle;

Diam_2 = MASS.Diam_2;
Diam_1 = MASS.Diam_1;

x_fair = MASS.TR.Fair.fair_length;

if strcmp(MASS.Tank2.OX.FLAG, 'Cyl') && strcmp(MASS.Tank2.FU.FLAG, 'Cyl')
% Top FU, bottom OX
L_fwrd_skirt_s2 = ((1/3)*Diam_2) + MASS.Tank2.FU.H_dome;
L_cyl_tank_FU2 = MASS.Tank2.FU.H_cyl;
L_connection_2 = ((1/4)*Diam_2) + (MASS.Tank2.OX.H_dome + MASS.Tank2.FU.H_dome);
L_cyl_tank_OX2 = MASS.Tank2.OX.H_cyl;
L_dwrd_skirt_s2 = ((1/3)*Diam_2) + MASS.Tank2.OX.H_dome;

x_fwrd_skirt_s2 = x_fair + L_fwrd_skirt_s2;
x_cyl_22 = x_fwrd_skirt_s2 + L_cyl_tank_FU2;
x_connection_2 = x_cyl_22 + L_connection_2;
x_cyl_21 = x_connection_2 + L_cyl_tank_OX2;
x_dwrd_skirt_s2 = x_cyl_21 + L_dwrd_skirt_s2;

R_tank2_FU = MASS.Tank2.FU.R_cyl_fuel;
R_tank2_OX = MASS.Tank2.OX.R_cyl_lox;

elseif strcmp(MASS.Tank2.OX.FLAG, 'Sphere') && strcmp(MASS.Tank2.FU.FLAG, 'Sphere')

   L_fwrd_skirt_s2 = (MASS.Tank2.FU.R_sphere_fuel/2);
   L_rem_sphere_FU2 = MASS.Tank2.FU.R_sphere_fuel;
   L_connection_2 = ((1/4)*Diam_2) + ((MASS.Tank2.OX.R_sphere_lox/2) + (MASS.Tank2.FU.R_sphere_fuel/2));
   L_rem_sphere_OX2 = MASS.Tank2.OX.R_sphere_lox;
   L_dwrd_skirt_s2 = (MASS.Tank2.OX.R_sphere_lox/2) + ((1/3)*Diam_2);

x_fwrd_skirt_s2 = x_fair + L_fwrd_skirt_s2;
x_sph_22 = x_fwrd_skirt_s2 + L_rem_sphere_FU2;
x_connection_2 = x_sph_22 + L_connection_2;
x_sph_21 = x_connection_2 + L_rem_sphere_OX2;
x_dwrd_skirt_s2 = x_sph_21 + L_dwrd_skirt_s2;

R_tank2_FU = MASS.Tank2.FU.R_sphere_fuel;
R_tank2_OX = MASS.Tank2.OX.R_sphere_lox;

end

COM.Stk2.x_end_2 = x_dwrd_skirt_s2 + L_nozzle2;  % new end of LV, 2^nd stack end

N_parts_s2 = 10; % 9 + half interstage
L_interstage = Diam_1;%L_nozzle2*(1+0.1);
M_cables_s2_tot = MASS.M_cables_distributed*(x_dwrd_skirt_s2 + (L_interstage/2));
M_cables_s2 = M_cables_s2_tot/N_parts_s2;

%% C.O.M. positioning for each component

% x=0 at tip

COM.Stk2.Pay.x_pay = MASS.TR.Fair.fair_length/2;

COM.Stk2.Fair.x_fair = MASS.TR.Fair.fair_length*(2/3);

COM.Stk2.Frw_struct2.x_parachute_2 = x_fair + L_fwrd_skirt_s2/2;

if strcmp(MASS.Tank2.OX.FLAG, 'Cyl') && strcmp(MASS.Tank2.FU.FLAG, 'Cyl')

COM.Stk2.Prop2_FU.x_Mp2_FU0 = x_fwrd_skirt_s2+ L_cyl_tank_FU2/2;
COM.Stk2.Prop2_FU.x_Mp2_FU = COM.Stk2.Prop2_FU.x_Mp2_FU0 - (((abs(M_FU - MASS.Tank2.FU.M_fuel))/(MASS.Tank2.FU.M_fuel))*COM.Stk2.Prop2_FU.x_Mp2_FU0);
COM.Stk2.Tank2_FU.x_tank2_FU = COM.Stk2.Prop2_FU.x_Mp2_FU0;

COM.Stk2.Prop2_OX.x_Mp2_OX0 = x_connection_2 + L_cyl_tank_OX2/2;
COM.Stk2.Prop2_OX.x_Mp2_OX = COM.Stk2.Prop2_OX.x_Mp2_OX0 - ((abs((M_LOX - MASS.Tank2.OX.M_lox))/(MASS.Tank2.OX.M_lox))*COM.Stk2.Prop2_OX.x_Mp2_OX0);
COM.Stk2.Tank2_OX.x_tank2_OX = COM.Stk2.Prop2_OX.x_Mp2_OX0;

COM.Stk2.Pipes2.x_pipes2 = x_cyl_22 + Tank2.FU.H_dome; % about interface between 2 tanks

COM.Stk2.Aft_Struct2.x_aft_struct2 = x_cyl_21 + L_dwrd_skirt_s2/2;

COM.Stk2.Engine2.x_engine2 = x_cyl_21 + ((L_dwrd_skirt_s2)*(2/3));

elseif strcmp( MASS.Tank2.OX.FLAG, 'Sphere') && strcmp( MASS.Tank2.FU.FLAG, 'Sphere')

COM.Stk2.Prop2_FU.x_Mp2_FU0 =  x_fwrd_skirt_s2+ L_rem_sphere_FU2/2;
COM.Stk2.Prop2_FU.x_Mp2_FU = COM.Stk2.Prop2_FU.x_Mp2_FU0 - ((abs((M_FU -MASS.Tank2.FU.M_fuel ))/( MASS.Tank2.FU.M_fuel))*COM.Stk2.Prop2_FU.x_Mp2_FU0);
COM.Stk2.Tank2_FU.x_tank2_FU = COM.Stk2.Prop2_FU.x_Mp2_FU0;

COM.Stk2.Prop2_OX.x_Mp2_OX0 = x_connection_2 + L_rem_sphere_OX2/2;
COM.Stk2.Prop2_OX.x_Mp2_OX = COM.Stk2.Prop2_OX.x_Mp2_OX0 - (((abs(M_LOX - MASS.Tank2.OX.M_lox))/( MASS.Tank2.OX.M_lox))*COM.Stk2.Prop2_OX.x_Mp2_OX0);
COM.Stk2.Tank2_OX.x_tank2_OX = COM.Stk2.Prop2_OX.x_Mp2_OX0;

COM.Stk2.Pipes2.x_pipes2 = x_sph_22 + MASS.Tank2.FU.R_sphere_fuel/2; % about interface between 2 tanks

COM.Stk2.Aft_Struct2.x_aft_struct2 = x_sph_21 + L_dwrd_skirt_s2/2;

COM.Stk2.Engine2.x_engine2 = x_sph_21 + ((L_dwrd_skirt_s2)*(2/3));

end

%% Mass association to every component of 2nd stack (as before, no change in their mass)

COM.Stk2.Pay.m_pay = 250;

COM.Stk2.Fair.m_fair =  MASS.TR.Fair.M_fair+ M_cables_s2;

COM.Stk2.Frw_struct2.m_parachute_2 = 0.07* MASS.TR.M.Ms2+ M_cables_s2 + MASS.M_avionics_s2/2+  MASS.M_struct2/2;

COM.Stk2.Tank2_FU.m_tank2_FU =  MASS.Tank2.FU.M_tot_tank_fuel+ M_cables_s2;

COM.Stk2.Tank2_OX.m_tank2_OX =  MASS.Tank2.OX.M_tot_tank_lox+ M_cables_s2;

COM.Stk2.Prop2_FU.m_prop2_FU0 =  MASS.Tank2.FU.M_fuel;
COM.Stk2.Prop2_FU.m_prop2_FU = M_FU;

COM.Stk2.Prop2_OX.m_prop2_OX0 =  MASS.Tank2.OX.M_lox;
COM.Stk2.Prop2_OX.m_prop2_OX = M_LOX;

COM.Stk2.Pipes2.m_pipes2 = 10 +  MASS.M_cables_distributed;

COM.Stk2.Aft_Struct2.m_aft_struct2 =  MASS.M_struct2/2+ M_cables_s2+  MASS.M_avionics_s2/2;

COM.Stk2.Engine2.m_engine2 =  MASS.Engine2.M+ M_cables_s2;

%% Computation of position of COM wrt nose and M_tb1 (at t = tb1)

% Vector assembly ( careful with associating right components)

COM_vec_m = [COM.Stk2.Pay.m_pay;COM.Stk2.Fair.m_fair;COM.Stk2.Frw_struct2.m_parachute_2;COM.Stk2.Prop2_FU.m_prop2_FU;COM.Stk2.Tank2_FU.m_tank2_FU;COM.Stk2.Prop2_OX.m_prop2_OX;COM.Stk2.Tank2_OX.m_tank2_OX;COM.Stk2.Pipes2.m_pipes2;COM.Stk2.Aft_Struct2.m_aft_struct2;COM.Stk2.Engine2.m_engine2];
COM_vec_x = [COM.Stk2.Pay.x_pay;COM.Stk2.Fair.x_fair;COM.Stk2.Frw_struct2.x_parachute_2;COM.Stk2.Prop2_FU.x_Mp2_FU;COM.Stk2.Tank2_FU.x_tank2_FU;COM.Stk2.Prop2_OX.x_Mp2_OX;COM.Stk2.Tank2_OX.x_tank2_OX;COM.Stk2.Pipes2.x_pipes2;COM.Stk2.Aft_Struct2.x_aft_struct2;COM.Stk2.Engine2.x_engine2];

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

X_COM = COM.Stk2.TR.x_LV;

%% Moments of Inertia

% First the distance wrt to the C.O.M. is computed: (tb1: at t=tb1, [-]: at t)
% [m]

MoI.Stk2.Pay.x_pay = abs(COM.Stk2.TR.x_LV - COM.Stk2.Pay.x_pay);

MoI.Stk2.Fair.x_fair = abs(COM.Stk2.TR.x_LV - COM.Stk2.Fair.x_fair);

MoI.Stk2.Frw_struct2.x_parachute_2 = abs(COM.Stk2.TR.x_LV - COM.Stk2.Frw_struct2.x_parachute_2);

MoI.Stk2.Prop2_FU.x_Mp2_FU = abs(COM.Stk2.TR.x_LV - COM.Stk2.Prop2_FU.x_Mp2_FU);

MoI.Stk2.Tank2_FU.x_tank2_FU= abs(COM.Stk2.TR.x_LV - COM.Stk2.Tank2_FU.x_tank2_FU);


MoI.Stk2.Prop2_OX.x_Mp2_OX = abs(COM.Stk2.TR.x_LV - COM.Stk2.Prop2_OX.x_Mp2_OX);

MoI.Stk2.Tank2_OX.x_tank2_OX= abs(COM.Stk2.TR.x_LV - COM.Stk2.Tank2_OX.x_tank2_OX);

MoI.Stk2.Pipes2.x_pipes2 = abs(COM.Stk2.TR.x_LV - COM.Stk2.Pipes2.x_pipes2);


MoI.Stk2.Aft_Struct2.x_aft_struct2= abs(COM.Stk2.TR.x_LV - COM.Stk2.Aft_Struct2.x_aft_struct2);

MoI.Stk2.Engine2.x_engine2= abs(COM.Stk2.TR.x_LV - COM.Stk2.Engine2.x_engine2);

if COM.Stk2.Prop2_FU.m_prop2_FU<=0

MoI.Stk2.Prop2_FU.x_Mp2_FU = 0;

end

if COM.Stk2.Prop2_OX.m_prop2_OX<=0

MoI.Stk2.Prop2_OX.x_Mp2_OX =0;

end

% Inertia moments of components (J0) are computed only for the tanks, 
% interstage, fairing while all other quantities are treated as point masses:

% All in [kg/m^2]:

MoI.Stk2.Fair.J0r_fair = (3/10)*MASS.TR.Fair.M_fair*((Diam_2/2)^2);

MoI.Stk2.Fair.J0y_fair = (3/5)*MASS.TR.Fair.M_fair*((MASS.TR.Fair.fair_length)^2) + (3/20)*MASS.TR.Fair.M_fair*((Diam_2/2)^2);

if strcmp(MASS.Tank2.OX.FLAG, 'Cyl') && strcmp(MASS.Tank2.FU.FLAG, 'Cyl')

t_2O = (MASS.Tank2.OX.t_lox)/(MASS.Tank2.OX.R_cyl_lox + MASS.Tank2.OX.t_lox);
t_2F =(MASS.Tank2.FU.t_fuel)/(MASS.Tank2.FU.R_cyl_fuel + MASS.Tank2.FU.t_fuel);

r_2O = MASS.Tank2.OX.R_cyl_lox + MASS.Tank2.OX.t_lox;
r_1O =  MASS.Tank2.OX.R_cyl_lox;

r_2F = MASS.Tank2.FU.R_cyl_fuel + MASS.Tank2.FU.t_fuel;
r_1F =  MASS.Tank2.FU.R_cyl_fuel;


MoI.Stk2.Tank2_OX.J0r_tank2_OX = (MASS.Tank2.OX.M_tot_tank_lox*((MASS.Tank2.OX.R_cyl_lox)^2)*(1 - t_2O + ((t_2O^2)/2))) + ((4/5)*MASS.Tank2.OX.M_spherical_cap_lox*(MASS.Tank2.OX.R_cyl_lox));
MoI.Stk2.Tank2_FU.J0r_tank2_FU = (MASS.Tank2.FU.M_tot_tank_fuel*((MASS.Tank2.FU.R_cyl_fuel)^2)*(1 - t_2F + ((t_2F^2)/2))) + ((4/5)*MASS.Tank2.FU.M_spherical_cap_fuel*(MASS.Tank2.FU.R_cyl_fuel));

MoI.Stk2.Tank2_OX.J0y_tank2_OX = (1/12)*(MASS.Tank2.OX.M_tot_tank_lox*((3*((r_2O^2) +(r_1O^2))) + MASS.Tank2.OX.H_cyl));
MoI.Stk2.Tank2_FU.J0y_tank2_FU = (1/12)*(MASS.Tank2.FU.M_tot_tank_fuel*((3*((r_2F^2) +(r_1F^2))) + MASS.Tank2.FU.H_cyl));

elseif strcmp(MASS.Tank2.OX.FLAG, 'Sphere') && strcmp(MASS.Tank2.FU.FLAG, 'Sphere')

MoI.Stk2.Tank2_OX.J0r_tank2_OX = (2/5)*(MASS.Tank2.OX.M_tot_tank_lox)*((((MASS.Tank2.OX.R_sphere_lox+ MASS.Tank2.OX.t_lox)^5)-((MASS.Tank2.OX.R_sphere_lox)^5))/(((MASS.Tank2.OX.R_sphere_lox+ MASS.Tank2.OX.t_lox)^3)-((MASS.Tank2.OX.R_sphere_lox)^3))); % [kg*m^2]

MoI.Stk2.Tank2_FU.J0r_tank2_FU = (2/5)*(MASS.Tank2.FU.M_tot_tank_fuel)*((((MASS.Tank2.FU.R_sphere_fuel+ MASS.Tank2.FU.t_fuel)^5)-((MASS.Tank2.FU.R_sphere_fuel)^5))/(((MASS.Tank2.FU.R_sphere_fuel+ MASS.Tank2.FU.t_fuel)^3)-((MASS.Tank2.FU.R_sphere_fuel)^3)));% [kg*m^2]

MoI.Stk2.Tank2_OX.J0y_tank2_OX = MoI.Stk2.Tank2_OX.J0r_tank2_OX;
MoI.Stk2.Tank2_FU.J0y_tank2_FU = MoI.Stk2.Tank2_FU.J0r_tank2_FU;

end

% Vector assembly (careful with associating right components)

%MoI_vec_J0 = [0;MoI.Stk2.Fair.J0_fair;0;0;MoI.Stk2.Tank2_FU.J0_tank2_FU;0;MoI.Stk2.Tank2_OX.J0_tank2_OX;0;0;0];
MoI_vec_J0y = [0;MoI.Stk2.Fair.J0y_fair;0;0;MoI.Stk2.Tank2_FU.J0y_tank2_FU;0;MoI.Stk2.Tank2_OX.J0y_tank2_OX;0;0;0];
MoI_vec_x = [MoI.Stk2.Pay.x_pay;MoI.Stk2.Fair.x_fair;MoI.Stk2.Frw_struct2.x_parachute_2;MoI.Stk2.Prop2_FU.x_Mp2_FU;MoI.Stk2.Tank2_FU.x_tank2_FU;MoI.Stk2.Prop2_OX.x_Mp2_OX;MoI.Stk2.Tank2_OX.x_tank2_OX;MoI.Stk2.Pipes2.x_pipes2;MoI.Stk2.Aft_Struct2.x_aft_struct2;MoI.Stk2.Engine2.x_engine2];

MoI_vec_xquad = MoI_vec_x.*MoI_vec_x;

MoI_vec_mxquad = dot(COM_vec_m,MoI_vec_xquad);

% Moments of inertia of our LL at t = tb1 + t2 and at t = tb1

MoI.Stk2.TR.Jy = MoI_vec_mxquad + sum(MoI_vec_J0y); % [kg m^2] t = tb1 + t2
 

J_yaw = MoI.Stk2.TR.Jy;
J_pitch = J_yaw;

%% Roll:

MoI.Stk2.Pay.r_pay = Diam_2/4;

MoI.Stk2.Fair.r_fair = Diam_2/2;

MoI.Stk2.Frw_struct2.r_parachute_2 =Diam_2/2;

MoI.Stk2.Prop2_FU.r_Mp2_FU = R_tank2_FU;

MoI.Stk2.Tank2_FU.r_tank2_FU= R_tank2_FU;

MoI.Stk2.Prop2_OX.r_Mp2_OX = R_tank2_OX;

MoI.Stk2.Tank2_OX.r_tank2_OX= R_tank2_OX;

MoI.Stk2.Pipes2.r_pipes2 = (R_tank2_OX+ R_tank2_FU)/2;

MoI.Stk2.Aft_Struct2.r_aft_struct2= Diam_2/2;

MoI.Stk2.Engine2.r_engine2= Diam_2/2;

MoI_vec_J0r = [0;MoI.Stk2.Fair.J0r_fair;0;0;MoI.Stk2.Tank2_FU.J0r_tank2_FU;0;MoI.Stk2.Tank2_OX.J0r_tank2_OX;0;0;0];

MoI_vec_r = [MoI.Stk2.Pay.r_pay;MoI.Stk2.Fair.r_fair;MoI.Stk2.Frw_struct2.r_parachute_2;MoI.Stk2.Prop2_FU.r_Mp2_FU;MoI.Stk2.Tank2_FU.r_tank2_FU;MoI.Stk2.Prop2_OX.r_Mp2_OX;MoI.Stk2.Tank2_OX.r_tank2_OX;MoI.Stk2.Pipes2.r_pipes2;MoI.Stk2.Aft_Struct2.r_aft_struct2;MoI.Stk2.Engine2.r_engine2];

MoI_vec_rquad = MoI_vec_r.*MoI_vec_r;

MoI_vec_mrquad = dot(COM_vec_m,MoI_vec_rquad);

% Moments of inertia of our LV at t = t1 + t2 and at t = 0

MoI.Stk2.TR.Jr = MoI_vec_mrquad + sum(MoI_vec_J0r); % [kg m^2] t = t1 + t2

J_roll = MoI.Stk2.TR.Jr;

end % function