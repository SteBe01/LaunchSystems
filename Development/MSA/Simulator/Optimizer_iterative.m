clear, clc
close all
clear dyn
clear stage_Separation
clear orbit_revolution

addpath(genpath("Functions"))
addpath(genpath("Functions_internal"))
addpath(genpath("Functions_events"))

full_flight = 1;
MaxIter = 250;
iter = 1;

[stages, params, init] = loadMission();

[T, Y, idxStage] = run_simulator(stages, params, init, full_flight);
prev_min = 290e3;
prev_max = 315e3;

while abs(norm(Y(end, 1:2)) - params.Re - 400e3) > 1e1 && iter < MaxIter
    current_val(iter) = (prev_max + prev_min)/2;
    params.pitch.final_altitude = current_val(iter);
    [T, Y, idxStage] = run_simulator(stages, params, init, full_flight);

    if norm(Y(end, 1:2)) - params.Re > 400e3
        prev_max = current_val(iter);
    elseif norm(Y(end, 1:2)) - params.Re < 400e3
        prev_min = current_val(iter);
    else
        break;
    end

    disp("Value used: " + num2str(current_val(iter)));
    iter = iter+1;
end