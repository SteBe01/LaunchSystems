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
b = [stages.stg1.m_prop/3 stages.stg2.m_prop/3 0 0 deg2rad(60) -deg2rad(15)]';
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
options.UseParallel = true;

init_guess = [0.1*stages.stg1.m_prop 0.1*stages.stg2.m_prop deg2rad(40)];

%%

[x, fval, exitFlag, output,lambda,grad,hessian] = fmincon(@(X) objFun([X(1) X(2) X(3)], stages, params, init),init_guess,A,b,Aeq,beq,lb,ub,nonlcon,options);

stages.stg1.m_prop_final = x(1);
stages.stg2.m_prop_final = x(2);
params.pitch.first_angle = x(3);

params.dispStat = true;
[T, Y, idxStage, parout] = run_simulator(stages, params, init, 1);
plotData(T, Y, params, parout, idxStage);

%% Auxiliary functions

function obj = objFun(x, stages, params, init)
    clear dyn
    clear stage_Separation

    stages.stg1.m_prop_final = x(1);
    stages.stg2.m_prop_final = x(2);
    params.pitch.first_angle = x(3);

    [~, Y, ~] = run_simulator(stages, params, init, 1);

    obj = abs(norm(Y(end,1:2)) - (400e3+params.Re))/1e3;
end

function [c, ceq] = nonlinconstr(x, stages, params, init)
    stages.stg1.m_prop_final = x(1);
    stages.stg2.m_prop_final = x(2);
    params.pitch.first_angle = x(3);

    [~, Y] = run_simulator(stages, params, init, 1);

    v_orbit = sqrt(398600/(400+params.Re*1e-3));

    c = [];
    ceq(1) = (Y(end,3)/1e3 - v_orbit);%/(v_orbit*1e3);
end














