clc; clearvars; close all
clear dyn
clear stage_Separation

        % % angle definition
        % if h < params.h1
        %     angle = params.angle1;
        % elseif h < params.h2
        %     angle = params.angle2;
        % elseif h < params.h3
        %     angle = params.angle3;
        % else
        %     angle = params.angle4;
        % end


addpath("Functions\");
addpath("Functions_events\");

[stages, params, init] = loadMission();
params.dispStat = false;

A = [1 -1 0 0 0 0 0 0
    0 1 -1 0 0 0 0 0];
b = [0 0]';
Aeq = [];
beq = [];
lb = [params.Re+11e3 params.Re+11e3 params.Re+11e3 0 0 0 0 0 0];
ub = [params.Re+450e3 params.Re+450e3 params.Re+450e3 deg2rad(60) deg2rad(60) deg2rad(60) stages.stg1.m_prop/7 stages.stg2.m_prop/7];
nonlcon = @(x) nonlinconstr(x, stages, params, init);

options = optimoptions("fmincon");
options.Display = "iter";
options.Algorithm = "active-set";
options.MaxIterations = 250;
options.MaxFunctionEvaluations = 5000;
options.FunctionTolerance = 1e-6;
options.ConstraintTolerance = 1e-6;

init_guess = [params.Re+100e3 params.Re+200e3 params.Re+300e3 deg2rad(45) deg2rad(30) deg2rad(10) stages.stg1.m_prop/10 stages.stg2.m_prop/10];
params.angle4 = 0;

%%

number_of_points = 1e3;
stpts = zeros(number_of_points, length(init_guess));

for ii = 1:number_of_points
    stpts(ii,1) = randi([params.Re+11e3 params.Re+399e3], 1,1);
    stpts(ii,2) = randi([stpts(ii,1) params.Re+399e3], 1,1);
    stpts(ii,3) = randi([stpts(ii,2) params.Re+399e3], 1,1);

    stpts(ii,6) = randi([deg2rad(0) 2000], 1,1);
    stpts(ii,5) = randi([stpts(ii,6) 4000], 1,1);
    stpts(ii,4) = randi([stpts(ii,5) 5000], 1,1);
    stpts(ii,6) = deg2rad(stpts(ii,6)*1e-2);
    stpts(ii,5) = deg2rad(stpts(ii,5)*1e-2);
    stpts(ii,4) = deg2rad(stpts(ii,4)*1e-2);

    stpts(ii,7) = randi([0 round(stages.stg1.m_prop*1e2/8)], 1,1)*1e-2;
    stpts(ii,8) = randi([0 round(stages.stg2.m_prop*1e2/8)], 1,1)*1e-2;
end




objFun1 = @(X) objFun([X(1) X(2) X(3) X(4) X(5) X(6) X(7) X(8)], stages, params, init);

opts = optimoptions(@fmincon, 'Algorithm', 'sqp');
newprob = createOptimProblem('fmincon', 'x0', init_guess, 'Aineq',A, 'bineq',b, 'lb',lb, 'ub',ub, 'objective',objFun1, 'nonlcon',nonlcon, 'options',opts);

gs = GlobalSearch;

useParallel = 1;

if useParallel
    ms = MultiStart(gs, 'UseParallel', true, 'Display', 'iter');
    pool = gcp('nocreate');
    if isempty(pool)
        parpool
    else
        warning("Using existing parpool")
    end
else
    ms = MultiStart(gs, 'UseParallel', false, 'Display', 'iter');
end

startpts = CustomStartPointSet(stpts);

[xcust, fcust] = run(ms, newprob, startpts);

%%

% [x, fval, exitFlag, output] = fmincon(@(X) objFun([X(1) X(2) X(3) X(4) X(5) X(6) X(7) X(8)], stages, params, init),init_guess,A,b,Aeq,beq,lb,ub,nonlcon,options);

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

    params.h1 = x(1);
    params.h2 = x(2);
    params.h3 = x(3);
    params.angle1 = x(4);
    params.angle2 = x(5);
    params.angle3 = x(6);
    stages.stg1.m_prop_final = x(7);
    stages.stg2.m_prop_final = x(8);

    [~, Y, idxStage] = run_simulator(stages, params, init, 1);

    obj = Y(idxStage,2);
end

function [c, ceq] = nonlinconstr(x, stages, params, init)
    stages.stg1.m_prop_final = x(1);
    stages.stg2.m_prop_final = x(2);
    params.pitch.first_angle = x(3);

    [~, Y] = run_simulator(stages, params, init, 1);

    v_orbit = sqrt(398600/(400+params.Re*1e-3));

    c(1) = Y(end,2) - 400e3;
    c(2) = Y(end,3) - v_orbit;
    ceq = [];
end

