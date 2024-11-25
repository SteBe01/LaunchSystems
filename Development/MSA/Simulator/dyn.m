function [dY, parout] = dyn(t,y, stage, params, current_stage)

    % Y(1) = x (horizontal position, inertial)
    % Y(2) = z (vertical position, inertial)
    % Y(3) = xDot (horizontal velocity, inertial)
    % Y(4) = zDot (vertical velocity, inertial)
    % Y(5) = theta (pitch angle, positive counterclockwise)
    % Y(6) = thetaDot (pitch angular velocity, positive counterclockwise)

    % Persistent data
    persistent turn_complete
    if isempty(turn_complete)
        turn_complete = false;
    end

    % Retrieve data from ode
    x = y(1);
    z = y(2);
    xDot = y(3);
    zDot = y(4);
    theta = y(5);
    thetaDot = y(6);

    velsNorm = norm([xDot zDot]);
    dcm = [cos(theta) -sin(theta); sin(theta) cos(theta)]';
    gamma = atan2(xDot, zDot);
    alpha = theta-gamma;

    if current_stage == 1
        t_wait = params.t_turn;
    else 
        t_wait = stage.t_ign;
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

    % % Compute AoA
    % vels_body = dcm*[xDot; zDot];
    % wind_body = dcm*params.wind_ned;
    % v_rel = vels_body - wind_body;
    % alpha = atan2(v_rel(2), v_rel(1));
    
    % Environment data
    g = params.g0/((1+z/Re)^2);                     % [m/s^2]   - Gravity acceleration taking into account altitude
    rho = getDensity(z);                            % [kg/m^3]  - Density at current altitude
    qdyn = 0.5*rho*velsNorm^2;                             % [Pa]      - Dynamic pressure
    S = pi*(stage.d^2/4);                           % [m^2]     - Rocket surface area
    
    % Aerodynamic forces
    D = qdyn*S*stage.Cd;                            % [N]       - Drag force acting on the rocket
    L = qdyn*S*stage.Cl;                            % [N]       - Lift force acting on the rocket
    
    % Thrust vectoring
    % if stage.useTVC
    %     t_turn = nan;
    %     if current_stage == 1 && t >= params.t_turn && abs(theta - params.theta_turn) > 1e-3 && ~turn_complete
    %         delta = stage.deltaMax;
    %     elseif current_stage == 1 && t >= params.t_turn && abs(theta - params.theta_turn) <= 1e-3 && ~turn_complete
    %         t_turn = t;
    %         turn_complete = true;
    %         % delta = -params.k1*theta - params.k2*thetaDot - params.k3*alpha;
    %         delta = 0;
    %     else
    %         % delta = -params.k1*theta - params.k2*thetaDot - params.k3*alpha;
    %         delta = 0;
    %     end
    % else
    %     delta = 0;
    % end

    if theta > 0
        delta = deg2rad(5);
    else
        delta = deg2rad(-5);
    end

    % Thrust
    if t > t_wait && t <= t_burn_tot + t_wait       % [N]       - Thrust force acting on the rocket
        T = stage.Thrust;
    else
        T = 0;
    end

    % Forces on the rocket in body frame
    F_i = -D*cos(alpha) - L*sin(alpha) + T*cos(delta) - m*g*cos(gamma);
    F_j = +D*sin(alpha) - L*cos(alpha) + T*sin(delta) + m*g*sin(gamma);

    % Moments on the rocket in body frame
    M_t = (D*sin(alpha)-L*cos(alpha))*(stage.xcp - stage.xcg) + T*sin(delta)*(stage.length - stage.xcg);

   % Forces in inertial frame
   Fxz = dcm'*[F_i; F_j];
   F_x = Fxz(1);
   F_z = Fxz(2);

    % Derivative vector
    dY = zeros(6, 1);

    dY(1) = xDot;
    dY(2) = zDot;
    dY(3) = F_x/m;
    dY(4) = F_z/m;
    dY(5) = thetaDot;
    dY(6) = M_t/stage.I;

    % Prepare output struct for ode recall
    if nargout > 1
        parout.qdyn = qdyn;
        parout.acc = reshape(Fxz, [1 2])/m;
        if exist("t_turn", 'var') && ~isnan(t_turn)
            parout.t_turn = t_turn;
        end
        parout.alpha = alpha;
    end
end

