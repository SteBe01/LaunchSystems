function [dY, parout] = dyn(t,y, stage, params, current_stage)

    % Y(1) = v (velocity, body frame)
    % Y(2) = gamma (trajectory path angle)
    % Y(3) = x (downrange)
    % Y(4) = h (altitude)
    % Y(5) = theta
    % Y(6) = thetaDot

    % Persistent data
    persistent turn_complete gamma_drop
    if isempty(turn_complete)
        turn_complete = false;
    end

    % Retrieve data from ode
    v = y(1);
    gamma = y(2);
    % x = y(3);
    h = y(4);
    % theta = y(5);
    thetaDot = y(6);

    if t < params.t_turn
        gamma_drop = gamma;
    end

    if current_stage == 1
        t_wait = params.t_turn;
    else 
        t_wait = stage.t_ign;
    end

    % Thrust vectoring
    if stage.useTVC
        t_turn = nan;
        if current_stage == 1 && t >= params.t_turn && (gamma - params.gamma_turn) < 1e-3 && ~turn_complete
            delta = stage.deltaMax;
        elseif current_stage == 1 && t >= params.t_turn && (gamma - params.gamma_turn) >= 1e-3 && ~turn_complete
            t_turn = t;
            turn_complete = true;
            delta = 0;
        else
            delta = 0;
        end
    else
        delta = 0;
    end

    % Retrieve data used multiple times 
    t_burn_tot = stage.t_burn_tot;
    Re = params.Re;

    if t <= t_wait
        m = stage.m0;
    elseif t > t_wait && t <= t_burn_tot + t_wait
        m = stage.m0 - stage.m_dot * (t-t_wait);
    else
        m = stage.m0 - stage.m_dot * t_burn_tot;
    end
    
    % Forces
    g = params.g0/((1+h/Re)^2);                     % [m/s^2]   - Gravity acceleration taking into account altitude
    rho = getDensity(h);                            % [kg/m^3]  - Density at current altitude
    qdyn = 0.5*rho*v^2;                             % [Pa]      - Dynamic pressure
    S = pi*(stage.d^2/4);                           % [m^2]     - Rocket surface area
    D = qdyn*S*stage.Cd;                            % [N]       - Drag force acting on the rocket
    L = qdyn*S*stage.Cl;                            % [N]       - Lift force acting on the rocket
    if t > t_wait && t <= t_burn_tot + t_wait       % [N]       - Thrust force acting on the rocket
        T = stage.Thrust;
    else
        T = 0;
    end

    % Normal acceleration
    acc_N = (T*sin(delta) + L - m*g*cos(gamma))/m;

    % Moments
    Ma = 0;
    Md = T*sin(delta)*1.5/stage.I;

    % Derivative vector
    dY = zeros(4, 1);

    if current_stage == 1 && t < t_wait
        dY(1) = -D/m;
        dY(3) = Re/(Re+h) * v;
        dY(4) = -g*t;
    else
        dY(1) = T/m*cos(delta) - D/m - g*sin(gamma);
        dY(2) = v*cos(gamma)/(Re+h) - g*cos(gamma)/v + T*sin(delta)/(m*v) + L/(m*v);
        dY(3) = Re/(Re+h) * v * cos(gamma);
        dY(4) = v * sin(gamma);
    end

    dY(5) = thetaDot;
    dY(6) = Md*1;

    % Pitch-up maneuver
    % if current_stage == 1 && (t >= params.t_turn && t <= params.t_turn+params.turn_duration)
    %     dY(2) = (params.gamma_turn-gamma_drop)*(pi/(2*params.turn_duration)*sin(pi*(t-params.t_turn)/params.turn_duration));
    % end

    % Prepare output struct for ode recall
    if nargout > 1
        parout.qdyn = qdyn;
        parout.gamma_dot = dY(2);
        parout.acc = [dY(1); acc_N];
        if exist("t_turn", 'var') && ~isnan(t_turn)
            parout.t_turn = t_turn;
        end
    end
end

