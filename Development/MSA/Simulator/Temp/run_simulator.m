function [T, Y, idxStage, parout] = run_simulator(stages, params, init)

    clear dyn stage_Separation

    % Trajectory propagation
    y0_stg1 = [init.x0 init.z0 init.vx0 init.vz0 init.theta0 init.thetaDot0 init.m_prop];

    options_stg1 = odeset('RelTol',1e-8, 'MaxStep', 0.1, 'Events', @(t, y) stage_Separation(t, y, stages.stg1));
    options_stg2 = odeset('RelTol', 1e-8, 'MaxStep', 0.1, 'Events', @(t, y) orbit_insertion(t, y));

    t_max = 1e4;
    
    %% First stage simulation
    [T1, Y1] = ode113(@(t,y) dyn(t, y, stages.stg1, params, 1), [0 100], y0_stg1, options_stg1);
    clear dyn

    %% Second stage simulation
    [T2, Y2] = ode113(@(t,y) dyn(t, y, stages.stg2, params, 2), [0 1], [Y1(end,1:end-1) stages.stg2.m_prop], options_stg2);
    clear dyn

    %% Retrieve data from ode

    T = [T1; T2+T1(end)];
    Y = [Y1; Y2];

    idxStage = length(T1);

    if nargout > 3

        parout_stg1 = recallOdeFcn(T1, Y1, stages.stg1, params, 1);
        parout_stg2 = recallOdeFcn(T2, Y2, stages.stg2, params, 2);

        parout.qdyn = [parout_stg1.qdyn; parout_stg2.qdyn];
        parout.acc = [parout_stg1.acc; parout_stg2.acc];
        parout.alpha = [parout_stg1.alpha; parout_stg2.alpha];
        parout.moment = [parout_stg1.moment; parout_stg2.moment];
        parout.dv_drag_vec = [parout_stg1.dv_drag_vec; parout_stg2.dv_drag_vec];
        parout.dv_grav_vec = [parout_stg1.dv_grav_vec; parout_stg2.dv_grav_vec];
        parout.delta_vec = [parout_stg1.delta; parout_stg2.delta];

        dv_drag_s1 = cumtrapz(parout_stg1.dv_drag_vec, T1);
        dv_grav_s1 = cumtrapz(parout_stg1.dv_grav_vec, T1);
        dv_thrust_s1 = stages.stg1.Isp*9.81*log(stages.stg1.m0/parout_stg1.m(end));
        parout.dv_s1 = dv_drag_s1(end) + dv_grav_s1(end) + dv_thrust_s1;

        dv_drag_s2 = cumtrapz(parout_stg2.dv_drag_vec, T2);
        dv_grav_s2 = cumtrapz(parout_stg2.dv_grav_vec, T2);
        dv_thrust_s2 = stages.stg2.Isp*9.81*log(stages.stg2.m0/parout_stg2.m(end));
        parout.dv_s2 = dv_drag_s2(end) + dv_grav_s2(end) + dv_thrust_s2;

        parout.coeffs = [parout_stg1.coeffs; parout_stg2.coeffs];

        parout.dcm = [parout_stg1.dcm];
        parout.F_in = [parout_stg1.F_in];

        parout.F_L_in = parout_stg1.F_L_in;
        parout.F_D_in = parout_stg1.F_D_in;
    end

end

