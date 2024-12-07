function [T, Y, idxStage, parout] = run_simulator_fs_reentry(stages, params, init)


    %% Final conditions of first stage

    clear dyn stage_Separation orbit_revolution
    t_max = 1e3;
    
    y0_stg1 = [init.x0 init.z0 init.vx0 init.vz0 init.theta0 init.thetaDot0 init.m_prop];

    options_stg1 = odeset('RelTol',1e-8, 'MaxStep', 0.1, 'Events', @(t, y) stage_Separation(t, y, stages.stg1, params));
    [~, Y1] = ode113(@(t,y) dyn(t, y, stages.stg1, params, 1), [0 t_max], y0_stg1, options_stg1);
    clear dyn


    %% First stage simulation

    clear dyn_fs_reentry
    t_max = 1e3;
        
    y0_stg1 = Y1(end,:);
    y0_stg1 = [y0_stg1 300];

    options_stg1 = odeset('RelTol',1e-8, 'MaxStep', 0.1, 'Events', @(t, y) touchdown(t, y, params));

    [T1, Y1, ~, ~, ie] = ode113(@(t,y) dyn_fs_reentry(t,y, stages, params), [0 t_max], y0_stg1, options_stg1);
    
    clear dyn

    T = T1;
    Y = Y1;

    
    %% Retrieve data from ode

    idxStage = length(T1);

    if nargout > 3

        parout_stg1 = recallOdeFcn(T1, Y1, stages.stg1, params, 1);

        parout.qdyn = parout_stg1.qdyn;
        parout.acc = parout_stg1.acc;
        parout.alpha = parout_stg1.alpha;
        parout.moment = parout_stg1.moment;
        parout.dv_drag_vec = parout_stg1.dv_drag_vec;
        parout.dv_grav_vec = parout_stg1.dv_grav_vec;
        parout.delta_vec = parout_stg1.delta;
        parout.m_vec = parout_stg1.m;

        dv_drag_s1 = cumtrapz(parout_stg1.dv_drag_vec, T1);
        dv_grav_s1 = cumtrapz(parout_stg1.dv_grav_vec, T1);
        dv_thrust_s1 = stages.stg1.Isp*9.81*log(stages.stg1.m0/parout_stg1.m(end));
        parout.dv_s1 = dv_drag_s1(end) + dv_grav_s1(end) + dv_thrust_s1;

        parout.coeffs = parout_stg1.coeffs;

        for ii = 1:idxStage
            parout.dcm(:,:,ii) = parout_stg1.dcm(:,:,ii);
        end
        parout.F_in = parout_stg1.F_in;

        parout.F_L_in = parout_stg1.F_L_in;
        parout.F_D_in = parout_stg1.F_D_in;
        parout.Thrust = parout_stg1.Thrust;

        parout.idxStg1 = length(parout_stg1.qdyn);
    end
    if ie == 2
        parout.dv_s1 = parout.dv_s1/T(end)*100*1e3;
    end
end

