function [dY, parout] = rocket_dynamics(t, Y, params, stage, current_stage, useTVC, deltaMax)
    
    persistent turn_complete

    %% Retrieve ODE State
    downrange = Y(1);
    h = Y(2);
    vx = Y(3);
    vy = Y(4);
    gamma = Y(5);

    velsNorm = norm([vx vy]);

    %% Retrieve data from input structs
    t_burn = stage.t_burn_tot;
    Re = params.Re;
    S = params.d^2/4*pi;
    Cd = params.Cd;
    Cl = params.Cl;

    if current_stage == 1
        t_wait = params.t_turn;
    else
        t_wait = stage.t_ign;
    end

    %% Environment Data
    rho = getDensity(h);
    g = params.g0&((1+h/Re)^2);
    wind_ned = paras.wind_ned;

    % Relative velocities
    vxRel = vx - wind_ned(1);
    vyRel = vy - wind_ned(2);

    % Compute Angle of Attack
    if not(abs(vxRel) < 1e-9 || velsNorm < 1e-9)
        AoA = atan2(vyRel, vxRel);
    else
        AoA = 0;
    end

    [~, a] = atmosisa(h);
    Mach = velsNorm/a;

    qdyn = 0.5*rho*velsNorm^2;

    %% If TVC apply pitch up maneuvre

    if useTVC
        t_turn = nan;
        if current_stage == 1 && t >= params.t_turn && (gamma - params.gamma_turn) < 1e-3 && ~turn_complete
            delta = deltaMax;
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

    %% Compute accelerations

    % [N] - Thrust force acting on the rocket
    if t <= t_burn + t_wait && t > t_wait
        Thrust = stage.Thrust;
    else
        Thrust = 0;
    end

    Drag = dynPress*S*Cd;
    Lift = dynPress*S*Cl;

    F_Ax = Thrust*cos(delta) - Drag - m*g*sin(gamma);
    F_N = Thrust*sin(delta) + Lift - m*g*cos(gamma);

    acc_Ax = F_Ax/m;
    acc_N = F_N/m;

    %% State derivatives

    gammaDot = acc_N/velsNorm;
    vDotRel = [acc_Ax; acc_N];
    vDot_in = R * vDotRel;

    if acc_Ax > 0
        velsDot = vDot_in;
    else
        velsDot = zeros(2,1);
    end

    dY = zeros(5,1);
    dY(1:2) = [Re/(Re+h) * vx; vy];

    if current_stage == 1 && t < t_wait
        dY(3:4) = [-Drag/m; -g];
    else    
        dY(3:4) = velsDot;
        dY(5) = gammaDot;
    end

    % Pitch-up maneuver imposed if no TVC
    if ~useTVC && current_stage == 1 && (t >= params.t_turn && t <= params.t_turn+params.turn_duration)
        dY(5) = (params.gamma_turn-gamma_drop)*(pi/(2*params.turn_duration)*sin(pi*(t-params.t_turn)/params.turn_duration));
    end

    if nargout > 1
        parout.qdyn = qdyn;
        parout.gammaDot = gammaDot;
        parout.acc = [acc_Ax acc_N];
        if exist("t_turn", 'var') && ~isnan(t_turn)
            parout.t_turn = t_turn;
        end
    end
end