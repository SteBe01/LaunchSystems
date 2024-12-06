%% Simulator

clear, clc
close all
clear dyn
clear stage_Separation
clear orbit_revolution

full_flight = 1;

addpath(genpath("Functions"))
addpath(genpath("Functions_events"))

[stages, params, init] = loadMission();

[T, Y, idxStage, parout] = run_simulator(stages, params, init, full_flight);
plotData(T, Y, params, parout, idxStage);

