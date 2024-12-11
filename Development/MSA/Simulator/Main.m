%% Simulator

clear clc
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

figure
subplot(2,1,1), hold on, grid on, title("$\dot Q$ over time", Interpreter="latex"), xlabel("Time [s]"), ylabel("$\dot Q$", Interpreter="latex")
plot(T, parout.Q)
subplot(2,1,2), hold on, grid on, title("Temperature over time"), xlabel("Time [s]"), ylabel("Temperature [K]")
plot(T, Y(:,8))


%% Second stage - reentry

clear, clc
close all
clear dyn_ss_reentry

addpath(genpath("Functions"))
addpath(genpath("Functions_internal"))
addpath(genpath("Functions_events"))

[stages, params, init] = loadMission();

[T, Y, idxStage, parout] = run_simulator_ss_reentry(stages, params);
plotData(T, Y, params, parout, idxStage);

figure
subplot(2,1,1), hold on, grid on, title("$\dot Q$ over time", Interpreter="latex"), xlabel("Time [s]"), ylabel("$\dot Q$", Interpreter="latex")
plot(T, parout.Q)
subplot(2,1,2), hold on, grid on, title("Temperature over time"), xlabel("Time [s]"), ylabel("Temperature [K]")
plot(T, Y(:,8))

