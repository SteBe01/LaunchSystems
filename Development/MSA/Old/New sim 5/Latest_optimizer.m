%%

clc; clearvars; close all
clear dyn
clear stage_Separation

addpath("Functions\");
addpath("Functions_events\");

[stages, params, init] = loadMission();
params.dispStat = false;

A = [];
b = []';
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = @(x) nonlinconstr(x, stages, params, init);

options = optimoptions("fmincon");
options.Display = "iter";
options.Algorithm = "active-set";
options.MaxIterations = 250;
options.MaxFunctionEvaluations = 5000;
options.FunctionTolerance = 1e-6;
options.ConstraintTolerance = 1e-6;

init_guess = [];

%%

[x, fval, exitFlag, output] = fmincon(@(X) objFun([X(1) X(2) X(3) X(4)], stages, params, init),init_guess,A,b,Aeq,beq,lb,ub,nonlcon,options);

%%

% stages.stg1.m_prop_final = x(1);
% stages.stg2.m_prop_final = x(2);
% params.pitch.first_angle = x(3);
% 
% params.dispStat = true;
% [T, Y, idxStage, parout] = run_simulator(stages, params, init);
% plotData(T, Y, params, parout, idxStage);

%% Auxiliary functions

function obj = objFun(x, stages, params, init)
    clear dyn
    clear stage_Separation

    stages.stg1.m_prop_final = x(1)*stages.stg1.m_prop;
    stages.stg2.m_prop_final = x(2)*stages.stg2.m_prop;
    
    params.pitch.first_angle = deg2rad(x(3));
    params.pitch.order = 2;
    params.pitch.initial_altitude = 11e3;
    params.pitch.final_altitude = x(4);

    [~, Y, idxStage] = run_simulator(stages, params, init, 1);

    obj = norm([Y(idxStage,2) Y(idxStage, 3)]);
end

function [c, ceq] = nonlinconstr(x, stages, params, init)
    stages.stg1.m_prop_final = x(1);
    stages.stg2.m_prop_final = x(2);
    params.pitch.first_angle = x(3);

    [~, Y] = run_simulator(stages, params, init, 1);

    v_orbit = sqrt(398600/(400+params.Re*1e-3));

    c(1) = abs(Y(end,2) - 400e3) - 100;
    c(2) = abs(Y(end,3) - v_orbit) - 100;
    ceq = [];
end














