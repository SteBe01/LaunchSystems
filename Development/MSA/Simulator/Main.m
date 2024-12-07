%% Simulator

clear, clc
close all
clear dyn
clear stage_Separation
clear orbit_revolution

addpath(genpath("Functions"))
addpath(genpath("Functions_internal"))
addpath(genpath("Functions_events"))

full_flight = 1;

[stages, params, init] = loadMission();

[T, Y, idxStage, parout] = run_simulator(stages, params, init, full_flight);
plotData(T, Y, params, parout, idxStage);


%% First stage - reentry

clear, clc
close all
clear dyn_fs_reentry

addpath(genpath("Functions"))
addpath(genpath("Functions_internal"))
addpath(genpath("Functions_events"))

[stages, params, init] = loadMission();

[T, Y, idxStage, parout] = run_simulator_fs_reentry(stages, params, init);
plotData(T, Y, params, parout, idxStage);

