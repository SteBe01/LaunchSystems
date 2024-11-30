clc; clearvars; close all
clear dyn
clear stage_Separation

addpath("Functions\");
addpath("Functions_events\");

[stages, params, init] = loadMission();
params.dispStat = false;

A = [1  0  0
     0  1  0
    -1  0  0
     0 -1  0
     0  0  1
     0  0 -1];
b = [stages.stg1.m_prop/3 stages.stg2.m_prop/3 0 0 pi/2 0]';
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

init_guess = [0.05*stages.stg1.m_prop 0.05*stages.stg2.m_prop deg2rad(60)];

%%

[x, fval, exitFlag, output] = fmincon(@(X) objFun([X(1) X(2) X(3)], stages, params, init),init_guess,A,b,Aeq,beq,lb,ub,nonlcon,options);

stages.stg1.m_prop_final = x(1);
stages.stg2.m_prop_final = x(2);
params.pitch.first_angle = x(3);

params.dispStat = true;
[T, Y, idxStage, parout] = run_simulator(stages, params, init);
plotData(T, Y, params, parout, idxStage);

%% Auxiliary functions

function obj = objFun(x, stages, params, init)
    clear dyn
    clear stage_Separation

    stages.stg1.m_prop_final = x(1);
    stages.stg2.m_prop_final = x(2);
    params.pitch.first_angle = x(3);

    [~, Y, idxStage] = run_simulator(stages, params, init);

    obj = norm([Y(idxStage,2) Y(idxStage, 3)]);
end

function [c, ceq] = nonlinconstr(x, stages, params, init)
    stages.stg1.m_prop_final = x(1);
    stages.stg2.m_prop_final = x(2);
    params.pitch.first_angle = x(3);

    [~, Y] = run_simulator(stages, params, init);

    c(1) = abs(Y(end,2) - 400e3) - 100;
    c(2) = abs(Y(end,3) - 7.6e3) - 100;
    ceq = [];
end














