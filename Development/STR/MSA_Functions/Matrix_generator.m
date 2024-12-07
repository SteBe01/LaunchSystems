clc; clearvars; close all

addpath("..\");

cd("..\..\MSA\Simulator\");
[stages, params, init] = loadMission();
cd("..\..\STR\MSA_Functions\")

h_it=load("h_it.mat");
M_it=load("M_it.mat");
M1=load("M1.mat");
M2=load("M2.mat");

n_choice = 75; % choice of case-> we chose 34

h_it_case = h_it.h_it(1,n_choice);
M_it_case=M_it.M_it(n_choice);

GEOMETRY = GEO(h_it_case,M_it_case,M1.M1,M2.M2);

N = 100;
outMat = zeros(N, 5, 2);

%% Stage 1
t_vec = linspace(0, stages.stg1.t_burn_tot, N)';

m_prop_left = zeros(N, 1);
xcg_vec = zeros(N, 1);
J_vec = zeros(N,3);

tic
for ii = 1:N
    m_prop_left(ii) = stages.stg1.m_prop - stages.stg1.m_dot(end)*stages.stg1.N_mot*t_vec(ii);
    [X_COM_t,Jyaw_t,Jpitch_t,Jroll_t] = MASS_PROPERTIES_1(m_prop_left(ii),GEOMETRY);
    xcg_vec(ii) = X_COM_t;
    J_vec(ii,:) = [Jyaw_t Jpitch_t Jroll_t];    
end
toc

outMat(:, :, 1) = [m_prop_left xcg_vec J_vec];

%% Stage 2
t_vec = linspace(0, stages.stg2.t_burn_tot, N)';

m_prop_left = zeros(N, 1);
xcg_vec = zeros(N, 1);
J_vec = zeros(N,3);

tic
for ii = 1:N
    m_prop_left(ii) = stages.stg2.m_prop - stages.stg2.m_dot(end)*stages.stg2.N_mot*t_vec(ii);
    [X_COM_t,Jyaw_t,Jpitch_t,Jroll_t] = MASS_PROPERTIES_2(m_prop_left(ii),GEOMETRY);
    xcg_vec(ii) = X_COM_t;
    J_vec(ii,:) = [Jyaw_t Jpitch_t Jroll_t];
end
toc

outMat(:, :, 2) = [m_prop_left xcg_vec J_vec];

%% Save output matrix

save("STR_mat.mat", "outMat");