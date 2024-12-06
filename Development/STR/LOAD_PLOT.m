clear;clc;close all;
%% RUNNARE PRIMA TRAIETTORIA E SCRIPT PIETRO:

stages=load('stages.mat');
stages = stages.stages;

% q = load("qdyn");
%q=max(q);
q=50*10^3;
FORCES.q = max(q);
% t_max_Q = find(Max_Q==qdyn);
% t_max_Q = t_max_Q(1,1);


FORCES.Cl_nose = stages.stg1.Cl;
FORCES.Cl_fin = stages.stg1.Cl;
FORCES.Cd = stages.stg1.Cd;

FORCES.T =stages.stg1.Thrust;
%FORCES.delta = delta_vec(t_max_Q);
FORCES.delta =0;
FORCES.alpha = 0;

FORCES.nx=7;
FORCES.nz=1.1;

% Script di Pietro: mass_opt_for_variable_dv.m
h_it=load("h_it.mat");
M_it=load("M_it.mat");
M1=load("M1.mat");
M2=load("M2.mat");

n_choice = 34; % choice of case-> we chose 34
h_it_case = h_it.h_it(1,n_choice);
M_it_case=M_it.M_it(n_choice);

GEOMETRY = GEO(h_it_case,M_it_case,M1.M1,M2.M2);

%% PLOT:
[STRUCT]=FINAL_STR_ANALYSIS(GEOMETRY,FORCES,2);
%[STRUCT]=FINAL_STR_ANALYSIS_T(GEOMETRY,FORCES,2);
