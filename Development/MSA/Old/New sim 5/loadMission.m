function [stages, params, init] = loadMission()

%% Display statistics data
params.dispStat = true;

%% Load MAT files

STR_mat = load("..\MAT_Files\STR_mat.mat").outMat;
stages.stg1.STR_mat = STR_mat(:,:,1);
stages.stg2.STR_mat = STR_mat(:,:,2);

PRP_mat = load("..\MAT_Files\PRP_mat.mat").throttle;
stages.stg1.Thrust = PRP_mat.TT;
stages.stg2.Thrust = PRP_mat.TT.*1.05;
stages.stg1.m_dot = PRP_mat.m_dot;
stages.stg2.m_dot = PRP_mat.m_dot;
stages.stg1.throttling = PRP_mat.manetta./100;
stages.stg2.throttling = PRP_mat.manetta./100;
stages.stg1.Pe = PRP_mat.p_e;
stages.stg2.Pe = PRP_mat.p_e;

wingSurf = 1;
AER_mat_path = strcat("../../AER/MSA/RISULTATI v2/", "S", num2str(wingSurf),"m^2/");
CD_mat = load(strcat(AER_mat_path, "data_CD_GEO2_S",num2str(wingSurf),"m^2.mat")).CD_mat;
CL_mat = load(strcat(AER_mat_path, "data_CL_GEO2_S",num2str(wingSurf),"m^2.mat")).CL_mat;
Xcp_mat = load(strcat(AER_mat_path, "data_Xcp_GEO2_S",num2str(wingSurf),"m^2.mat")).Xcp_mat;

AER_mat = cat(4, CD_mat, CL_mat, Xcp_mat);
AER_mat = permute(AER_mat, [4 1 2 3]);
params.coeffs = griddedInterpolant({(1:3)', (0.3:0.1:8)', (0:5:50)', (11000:500:84000)'}, AER_mat, "linear", "nearest");
params.CD_mult = 0.2;     % Multiplier for CD
params.CL_mult = 0.01;     % Multiplier for CL

%% General parameters
stages.stg1.MR = 3.4983;
stages.stg2.MR = 6.0115;
stages.stg1.Isp = 324.5;
stages.stg2.Isp = 343;
stages.stg1.A_eng = 0.0953;
stages.stg2.A_eng = 0.0953;
stages.stg1.N_mot = 8;
stages.stg2.N_mot = 1;
stages.stg1.m0 = 16.7e3;
stages.stg1.m_prop = stages.stg1.m0 * (1 - 1/stages.stg1.MR);
stages.stg1.t_burn_tot = stages.stg1.m_prop/(stages.stg1.m_dot(end)*stages.stg1.N_mot);
stages.stg2.m0 = stages.stg1.m0 - stages.stg1.m_prop - 1400;
stages.stg2.m_prop = stages.stg2.m0 * (1 - 1/stages.stg2.MR);
stages.stg2.t_burn_tot = stages.stg2.m_prop/(stages.stg2.m_dot(end)*stages.stg2.N_mot);

stages.stg1.m_prop_final = 0.0*stages.stg1.m_prop;
stages.stg2.m_prop_final = 0.05*stages.stg2.m_prop;

stages.stg1.d = 1.1636;
stages.stg2.d = 1.1636;
stages.stg1.length = 21.71;
stages.stg2.length = 8.16;

% Environment data
params.g0 = 9.81;
params.Re = 6378000;

% Initial conditions
init.x0 = 0;                                    % [m]
% init.z0 = 400e3+params.Re;                       % [m]
init.z0 = 11e3+params.Re;                       % [m]
init.vx0 = 200;                                 % [m/s]
% init.vx0 = sqrt(398600e9/(400e3+params.Re));                                 % [m/s]
init.vz0 = 0;                                   % [m/s]
init.theta0 = atan2(init.vz0,init.vx0);         % [rad]
init.thetaDot0 = deg2rad(0);                    % [rad/s]
init.m_prop = stages.stg1.m_prop;               % [kg]

% Stage 1 PID data
stages.stg1.k1 = deg2rad(3*3);
stages.stg1.k2 = 0.8;
stages.stg1.k3 = 0*3;
% Stage 2 PID Data
stages.stg2.k1 = deg2rad(9);
stages.stg2.k2 = 0.08;
stages.stg2.k3 = 0;

%% hardcoded data

% Controller frequency
stages.stg1.u_freq = 100;
stages.stg2.u_freq = 100;

% Pitch maneuver
params.t_turn = 1;                      % [s]       - Initial maneuver time
params.theta_turn = deg2rad(45);        % [rad]     - Initial flight path angle

% MECO to stage separation wait time
stages.stg1.t_wait = 5;

% Stage wait time before ignition
stages.stg2.t_ign = 3;

% TVC Parameters
stages.stg1.useTVC = true;
stages.stg2.useTVC = false;
stages.stg1.deltaMax = deg2rad(7);
stages.stg2.deltaMax = deg2rad(5);

%% pitch program

params.pitch.first_angle = deg2rad(40);
params.pitch.order = 1;
params.pitch.initial_altitude = 11e3;
params.pitch.final_altitude = 500e3;

end

