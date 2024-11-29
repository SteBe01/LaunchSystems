%% Simulator

clear, clc
close all
clear dyn
clear stage_Separation

addpath(genpath("Functions"))
addpath(genpath("Functions_events"))

[stages, params, init] = loadMission();

[T, Y, parout] = run_simulator_first_stage(stages, params, init);
plotData_first_stage(T, Y, params, parout);

