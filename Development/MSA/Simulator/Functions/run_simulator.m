function [T, Y, idxStage, parout] = run_simulator(stages, params, init, full_flight)

    clear dyn stage_Separation orbit_revolution

    T_orb = 2*pi*sqrt(6778^3/398600);
    t_max = 1e3;
    
    %% First stage simulation
    y0_stg1 = [init.x0 init.z0 init.vx0 init.vz0 init.theta0 init.thetaDot0 init.m_prop];

    if full_flight
        options_stg1 = odeset('RelTol',1e-8, 'MaxStep', 0.1, 'Events', @(t, y) stage_Separation(t, y, stages.stg1, params));
    else
        options_stg1 = odeset('RelTol',1e-8, 'MaxStep', 0.1, 'Events', @(t, y) touchdown(t, y, params));
        warning("Mass not updated (considering whole rocket, no stage separation)")
    end
    [T1, Y1, ~, ~, ie] = ode113(@(t,y) dyn(t, y, stages.stg1, params, 1), [0 t_max], y0_stg1, options_stg1);
    clear dyn

    T = T1;
    Y = Y1;

    %% Second stage simulation
    if full_flight && ie~=2
        options_stg2 = odeset('RelTol', 1e-8, 'MaxStep', 0.1);%, 'Events', @(t, y) orbit_insertion(t, y)); %@(t,y) orbit_revolution(t, y, params));
        [T2, Y2] = ode113(@(t,y) dyn(t, y, stages.stg2, params, 2), [0 t_max], [Y1(end,1:end-1) stages.stg2.m_prop], options_stg2);
        clear dyn

        T = [T; T2+T1(end)];
        Y = [Y; Y2];
    end

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

        if full_flight && ie~=2
            parout_stg2 = recallOdeFcn(T2, Y2, stages.stg2, params, 2, Y2(end, 7));
    
            parout.qdyn = [parout.qdyn; parout_stg2.qdyn];
            parout.acc = [parout.acc; parout_stg2.acc];
            parout.alpha = [parout.alpha; parout_stg2.alpha];
            parout.moment = [parout.moment; parout_stg2.moment];
            parout.dv_drag_vec = [parout.dv_drag_vec; parout_stg2.dv_drag_vec];
            parout.dv_grav_vec = [parout.dv_grav_vec; parout_stg2.dv_grav_vec];
            parout.delta_vec = [parout.delta_vec; parout_stg2.delta];
            parout.m_vec = [parout.m_vec; parout_stg2.m];
     
            dv_drag_s2 = cumtrapz(parout_stg2.dv_drag_vec, T2);
            dv_grav_s2 = cumtrapz(parout_stg2.dv_grav_vec, T2);
            dv_thrust_s2 = stages.stg2.Isp*9.81*log(stages.stg2.m0/parout_stg2.m(end));
            parout.dv_s2 = dv_drag_s2(end) + dv_grav_s2(end) + dv_thrust_s2;
    
            parout.coeffs = [parout.coeffs; parout_stg2.coeffs];
    
            for ii = 1:length(T2)
                parout.dcm(:,:,idxStage+ii) = parout_stg2.dcm(:,:,ii);
            end
            parout.F_in = [parout.F_in; parout_stg2.F_in];
    
            parout.F_L_in = [parout.F_L_in; parout_stg2.F_L_in];
            parout.F_D_in = [parout.F_D_in; parout_stg2.F_D_in];
            parout.Thrust = [parout.Thrust; parout_stg2.Thrust];
        end
        parout.idxStg1 = length(parout_stg1.qdyn);
    end
    if ie == 2
        parout.dv_s1 = parout.dv_s1/T(end)*100*1e3;
    end
end

