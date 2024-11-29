function [T, Y, parout] = run_simulator_first_stage(stages, params, init)

    % Trajectory propagation
    y0_stg1 = [init.x0 init.z0 init.vx0 init.vz0 init.theta0 init.thetaDot0 init.m_prop];

    options_stg1 = odeset('RelTol',1e-8, 'MaxStep', 0.1, 'Events', @(t, y) touchdown(t, y));
    t_max = 1e4;
    
    %% First stage simulation
    [T, Y] = ode113(@(t,y) dyn(t, y, stages.stg1, params, 1), [0 t_max], y0_stg1, options_stg1);
    clear dyn

    %% Retrieve data from ode

    parout_stg1 = recallOdeFcn(T, Y, stages.stg1, params, 1);

    parout.qdyn = [parout_stg1.qdyn];
    parout.acc = [parout_stg1.acc];
    parout.alpha = [parout_stg1.alpha];
    parout.moment = [parout_stg1.moment];
    parout.dv_drag_vec = [parout_stg1.dv_drag_vec];
    parout.dv_grav_vec = [parout_stg1.dv_grav_vec];
    parout.delta_vec = [parout_stg1.delta];

    dv_drag_s1 = cumtrapz(parout_stg1.dv_drag_vec, T);
    dv_grav_s1 = cumtrapz(parout_stg1.dv_grav_vec, T);
    dv_thrust_s1 = stages.stg1.Isp*9.81*log(stages.stg1.m0/parout_stg1.m(end));
    parout.dv_s1 = dv_drag_s1(end) + dv_grav_s1(end) + dv_thrust_s1;

end

%% Functions

function [value, isterminal, direction] = touchdown(~, y)
    value = y(2) - 25e3;
    isterminal = 1;
    direction = -1;
end

