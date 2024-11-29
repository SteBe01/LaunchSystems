%% Simulator

clear, clc
close all
clear dyn
clear stage_Separation

[stages, params, init] = loadMission();

[T, Y, parout, idxStage] = run_simulator(stages, params, init);
plotData(T, Y, params, parout, idxStage);

