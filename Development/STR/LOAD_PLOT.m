clear;clc;close all;
%% RUN TRAJECTORY AND PIETRO'S SCRIPT:

%% Trajectory:

% stages=load('stages.mat');
% stages = stages.stages;

% q = load("qdyn");
%q=max(q);
q=50*10^3;
FORCES.q = max(q);
% t_max_Q = find(Max_Q==qdyn);
% t_max_Q = t_max_Q(1,1);


% FORCES.Cl=stages.stg1.Cl;
% FORCES.Cd = stages.stg1.Cd;
FORCES.Cl=1.3;
FORCES.Cd = 1.3;

%FORCES.T =stages.stg1.Thrust;
%FORCES.delta = delta_vec(t_max_Q);
FORCES.delta =0;
FORCES.alpha = 0;

FORCES.nx=6;
FORCES.nz=1;

%% Script di Pietro: mass_opt_for_variable_dv.m
h_it=load("h_it.mat");
M_it=load("M_it.mat");

n_choice = 34; % choice of case-> we chose 34

h_it_case = h_it.h_it(1,n_choice);
M_it_case=M_it.M_it(n_choice);

GEOMETRY = GEO(h_it_case,M_it_case);

M_abs = abs(GEOMETRY.M01-GEOMETRY.m1-GEOMETRY.m2-GEOMETRY.m3-GEOMETRY.m4-GEOMETRY.m5-GEOMETRY.m6-GEOMETRY.m7-GEOMETRY.m8-GEOMETRY.m9-GEOMETRY.m10);

M_prop_vec_1 = linspace(GEOMETRY.m_prop_1,0,1000);
M_prop_vec_2 = linspace(GEOMETRY.m_prop_2,0,1000);

x_com1=zeros(length(M_prop_vec_1),1);
x_com2=zeros(length(M_prop_vec_1),1);
x_com2_nof=zeros(length(M_prop_vec_1),1);
x_com1_nof=zeros(length(M_prop_vec_1),1);
for i=1:length(M_prop_vec_1)

[x_com1(i)] = MASS_PROPERTIES_1(M_prop_vec_1(i),GEOMETRY);
[x_com2(i)] = MASS_PROPERTIES_2(M_prop_vec_2(i),GEOMETRY);
[x_com2_nof(i)] = MASS_PROPERTIES_2_NO_F(M_prop_vec_2(i),GEOMETRY);
[x_com1_nof(i)] = MASS_PROPERTIES_1_NO_F(M_prop_vec_1(i),GEOMETRY);

end

[X_COM10] = MASS_PROPERTIES_1(GEOMETRY.m_prop_1,GEOMETRY);
[X_COM10F] = MASS_PROPERTIES_1_NO_F(GEOMETRY.m_prop_1,GEOMETRY);
[X_COM20] = MASS_PROPERTIES_2(GEOMETRY.m_prop_2,GEOMETRY);
[X_COM20F] = MASS_PROPERTIES_2_NO_F(GEOMETRY.m_prop_1,GEOMETRY);

% plot(M_prop_vec_1,x_com1,'Color','b');
% hold on;
% plot(M_prop_vec_2,x_com2,'Color','r');
% hold on;
% plot(GEOMETRY.m_prop_1,X_COM10,'o');
% hold on;
% plot(GEOMETRY.m_prop_2,X_COM20,'o');
% plot(M_prop_vec_1,x_com1_nof,'Color','c');
% hold on;
% plot(M_prop_vec_2,x_com2_nof,'Color','#D95319');
% xlabel('Mass of Propellant [kg]');
% ylabel('X coord of COM wrt to nose [m]');
% legend('COM1','COM2','COM10','COM20','COM1 NO F','COM2 NO F','Location','eastoutside');
% grid on;
%M_prop_t1=GEOMETRY.m_prop_1;
%M_prop_t1=0;
%M_prop_t2=0;
%M_prop_t2 = GEOMETRY.m_prop_2;
%[X_COM_t1,Jyaw_t1,Jpitch_t1,Jroll_t1,GEOMETRY] = MASS_PROPERTIES_1(M_prop_t1,GEOMETRY);
%[X_COM_t2,Jyaw_t2,Jpitch_t2,Jroll_t2,GEOMETRY] = MASS_PROPERTIES_2(M_prop_t2,GEOMETRY);



% 
%% PLOT:
SaveFlag=1;
FORCES.T=M_it_case.stg1.Thrust;
%[STRUCT]=FINAL_STR_ANALYSIS(GEOMETRY,FORCES,2);
% % 
%[STRUCT]=FINAL_STR_ANALYSIS_T(GEOMETRY,FORCES,1,SaveFlag,h_it_case,M_it_case); 
%[STRUCT]=FINAL_STR_ANALYSIS_T(GEOMETRY,FORCES,2,SaveFlag,h_it_case,M_it_case); 
[STRUCT]=FINAL_STR_ANALYSIS_T_PLOT(GEOMETRY,FORCES,2,SaveFlag,h_it_case,M_it_case); 

 %% ATTACHMENTS:
PlotFlag =1;
LOAD.Cl=FORCES.Cl;
LOAD.Cd=FORCES.Cd;
LOAD.nx_c = 2;
LOAD.nz_c = 2;
%[CLAMP] = Attachments_FINAL(GEOMETRY,LOAD,PlotFlag,SaveFlag,h_it_case,M_it_case);
[CLAMP] = Attachments_FINAL_PLOT(GEOMETRY,LOAD,PlotFlag,SaveFlag,h_it_case,M_it_case);