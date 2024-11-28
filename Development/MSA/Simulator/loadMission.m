function [stages, params, init] = loadMission()

% General parameters
stages.stg1.MR = 5.225;
stages.stg2.MR = 4.022;
stages.stg1.Isp = 328;
stages.stg2.Isp = 343;
stages.stg1.Thrust = 27.4e3*9;
stages.stg2.Thrust = 31e3;
stages.stg1.m0 = 19.2e3;

stages.stg1.d = 1.8;
stages.stg2.d = 1.5;
stages.stg1.Cd = 0.8;
stages.stg2.Cd = 0.5;
stages.stg1.Cl = 2;
stages.stg2.Cl = 0.0;
stages.stg1.I = 3.579e5;
stages.stg2.I = 4.07e3;
stages.stg1.xcp = 16;
stages.stg2.xcp = 4;
stages.stg1.xcg = 7.5;
stages.stg2.xcg = 2; 
stages.stg1.length = 20.16;
stages.stg2.length = 4.6;

% Initial conditions
init.x0 = 0;                                    % [m]
init.z0 = 11e3;                                 % [m]
init.vx0 = 200;                                 % [m/s]
init.vz0 = 0;                                   % [m/s]
init.theta0 = atan2(init.vz0,init.vx0);         % [rad]
init.thetaDot0 = deg2rad(0);                    % [rad/s]

% Environment data
params.g0 = 9.81;
params.Re = 6378000;

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
params.t_turn = 5;                      % [s]       - Initial maneuver time
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

end

