function [stages, params] = loadMission()

% General parameters
stages.stg1.n = 5.53;
stages.stg2.n = 2.93;
stages.stg1.Isp = 280;
stages.stg2.Isp = 298;
stages.stg1.Thrust = 200e3;
stages.stg2.Thrust = 27.5e3;
stages.stg1.m0 = 13.5e3;

stages.stg1.d = 1.8;
stages.stg2.d = 1.5;
stages.stg1.Cd = 0.5;
stages.stg2.Cd = 0.5;
stages.stg1.Cl = 0.3;
stages.stg2.Cl = 0.0;
stages.stg1.I = 20e4;
stages.stg2.I = 5e3;

params.g0 = 9.81;
params.Re = 6378000;
params.h0 = 11.9e3;
% params.h0 = 400e3;
params.v0 = 200;
% params.v0 = 7.6686e3;
params.gamma0 = deg2rad(0);
params.turn_duration = 60;   % [s]
params.h_stage = 95e3;
params.wind_ned = [0; 0];

%% hardcoded data

% Pitch maneuver
params.t_turn = 5;                      % [s]       - Initial maneuver time
params.gamma_turn = deg2rad(45);        % [rad]     - Initial flight path angle

% MECO to stage separation wait time
stages.stg1.t_wait = 5;

% Stage wait time before ignition
stages.stg2.t_ign = 3;

% TVC Parameters
stages.stg1.useTVC = true;
stages.stg2.useTVC = false;
stages.stg1.deltaMax = deg2rad(7.5);
stages.stg2.deltaMax = 0;

end

