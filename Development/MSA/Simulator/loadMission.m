function [stages, params, init] = loadMission()

% General parameters
stages.stg1.n = 5.53;
stages.stg2.n = 2.93;
stages.stg1.Isp = 280;
stages.stg2.Isp = 298;
stages.stg1.Thrust = 27.4e3*9;
stages.stg2.Thrust = 27.5e3;
stages.stg1.m0 = 15e3;

stages.stg1.d = 1.8;
stages.stg2.d = 1.5;
stages.stg1.Cd = 0.5;
stages.stg2.Cd = 0.5;
stages.stg1.Cl = 0.0;
stages.stg2.Cl = 0.0;
stages.stg1.I = 4e5;
stages.stg2.I = 5e4;
stages.stg1.xcp = 10;
stages.stg2.xcp = 5;
stages.stg1.xcg = 6;
stages.stg2.xcg = 3; 
stages.stg1.length = 12;
stages.stg2.length = 6;

% Initial conditions
init.x0 = 0;                    % [m]
init.z0 = 11e3;               % [m]
init.vx0 = 200;                 % [m/s]
init.vz0 = 0;                   % [m/s]
init.theta0 = atan2(init.vz0,init.vx0);     % [rad]
init.thetaDot0 = deg2rad(0);    % [rad/s]

% Environment data
params.g0 = 9.81;
params.Re = 6378000;
params.turn_duration = 60;   % [s]
params.h_stage = 95e3;
params.wind_ned = [0; 0];

params.k1 = 0;
params.k2 = 0.8;
params.k3 = 3.614;

%% hardcoded data

% Pitch maneuver
params.t_turn = 5;                      % [s]       - Initial maneuver time
params.theta_turn = deg2rad(45);        % [rad]     - Initial flight path angle

% MECO to stage separation wait time
stages.stg1.t_wait = 5;

% Stage wait time before ignition
stages.stg2.t_ign = 3;

% TVC Parameters
stages.stg1.useTVC = true;
stages.stg2.useTVC = false;
stages.stg1.deltaMax = deg2rad(12);
stages.stg2.deltaMax = 0;

end

