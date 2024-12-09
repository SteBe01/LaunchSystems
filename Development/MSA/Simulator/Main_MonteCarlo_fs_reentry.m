%% Montecarlo simulations

clear, clc
close all
clear dyn
clear stage_Separation
clear orbit_revolution

rng default
N_sim = 100; % number of simulations

addpath(genpath("Functions"))
addpath(genpath("Functions_internal"))
addpath(genpath("Functions_events"))

full_flight = 1;

[nom_stages, nom_params, nom_init] = loadMission();
nom_params.dispStat = false;

%% Setup random data

mu_z0 = nom_init.z0;
sigma_z0 = 200/3; % 200 m at 3 sigma
z0_vec = normrnd(mu_z0, sigma_z0, N_sim, 1);

mu_vx0 = nom_init.vx0;
sigma_vx0 = 50/3; % 50 m/s at 3 sigma
vx0_vec = normrnd(mu_vx0, sigma_vx0, N_sim, 1);

mu_vz0 = nom_init.vz0;
sigma_vz0 = 50/3; % 50 m/s at 3 sigma
vz0_vec = normrnd(mu_vz0, sigma_vz0, N_sim, 1);

mu_t = 1;
sigma_t = 3/100;
thrust_percentage_stg1 = normrnd(mu_t,sigma_t,N_sim,1);

mu_mp = 1;
sigma_mp = 0.05/3;
prop_mass_percentage_stg1 = normrnd(mu_mp, sigma_mp, N_sim, 1);

mu_cl_mult = 1;
sigma_cl_mult = 0.1;
cl_mult = normrnd(mu_cl_mult, sigma_cl_mult, N_sim, 1);

mu_cd_mult = 1;
sigma_cd_mult = 0.1;
cd_mult = normrnd(mu_cd_mult, sigma_cd_mult, N_sim, 1);

mu_h1_offset = 0;
sigma_h1_offset = 1000/3;
h1_offset = normrnd(mu_h1_offset, sigma_h1_offset, N_sim, 1);

mu_h2_offset = 0;
sigma_h2_offset = 500/3;
h2_offset = normrnd(mu_h2_offset, sigma_h2_offset, N_sim, 1);

mu_h3_offset = 0;
sigma_h3_offset = 200/3;
h3_offset = normrnd(mu_h3_offset, sigma_h3_offset, N_sim, 1);

%% Prepare output vectors

saveSim = cell(N_sim, 1);
uncertVals = cell(N_sim, 1);

%% Run simualtions
tic
parfor ii = 1:N_sim
% for ii = 1:1
    % clear dyn
    % clear stage_Separation
    % clear orbit_revolution

    stages = nom_stages;
    params = nom_params;
    init = nom_init;

    init.z0 = z0_vec(ii);
    init.vx0 = vx0_vec(ii);
    init.vz0 = vz0_vec(ii);
    init.theta0 = atan2(init.vz0,init.vx0);
    stages.stg1.Thrust = stages.stg1.Thrust * thrust_percentage_stg1(ii);
    stages.stg1.m_prop = stages.stg1.m_prop * prop_mass_percentage_stg1(ii);
    stages.stg1.t_burn_tot = stages.stg1.m_prop/(stages.stg1.m_dot(end)*stages.stg1.N_mot);

    params.CL_mult = cl_mult(ii);
    params.CD_mult = cd_mult(ii);

    params.fs_h1 = nom_params.fs_h1 + h1_offset(ii);
    params.fs_h2 = nom_params.fs_h2 + h2_offset(ii);
    params.fs_h3 = nom_params.fs_h3 + h3_offset(ii);

    uncertVals{ii} = [z0_vec(ii); vx0_vec(ii); vz0_vec(ii); init.theta0;
                        thrust_percentage_stg1(ii); prop_mass_percentage_stg1(ii); 
                        cl_mult(ii); cd_mult(ii);
                        h1_offset(ii); h2_offset(ii); h3_offset(ii)];

    [T, Y, idxStage, parout] = run_simulator_fs_reentry(stages, params, init);
    saveSim{ii}.t = T;
    saveSim{ii}.Y = Y;
    saveSim{ii}.idxStage = idxStage;
    saveSim{ii}.parout = parout;

    fprintf("Completed simulation %d of %d\n", ii, N_sim);

end
totalTime = toc;
%%

fprintf("All simulations completed, total time required: %.3f minutes \n", totalTime/60);
% [T, Y, idxStage, parout] = run_simulator(stages, params, init, full_flight);
% plotData(T, Y, nom_params, parout, idxStage);