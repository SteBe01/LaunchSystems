function [stages, params] = loadMission()

% General parameters
stages.stg1.n = 5.53;
stages.stg2.n = 2.93;
stages.stg1.Isp = 280;
stages.stg2.Isp = 298;
stages.stg1.Thrust = 252e3;
stages.stg2.Thrust = 27.5e3;
stages.stg1.m0 = 15e3;

stages.stg1.d = 1.8;
stages.stg2.d = 1.5;
stages.stg1.Cd = 0.5;
stages.stg2.Cd = 0.5;
stages.stg1.Cl = 0.5;
stages.stg2.Cl = 0.0;

params.g0 = 9.81;
params.Re = 6378000;
params.h0 = 11.9e3;
params.v0 = 400;
params.gamma0 = deg2rad(0);
params.turn_duration = 60;   % [s]
params.h_stage = 95e3;

%% hardcoded data

% Pitch maneuver
params.t_turn = 5;                      % [s]       - Initial maneuver time
params.gamma_turn = deg2rad(30);        % [rad]     - Initial flight path angle

% MECO to stage separation wait time
stages.stg1.t_wait = 5;

% Stage wait time before ignition
stages.stg2.t_ign = 3;

end

