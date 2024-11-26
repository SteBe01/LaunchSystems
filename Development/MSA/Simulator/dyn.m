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
    rot_angle = theta - pi/2;
    dcm = [cos(rot_angle) -sin(rot_angle); sin(rot_angle) cos(rot_angle)];
    gamma = atan2(zDot, xDot);
    alpha = -theta+gamma;

    % if current_stage == 1
    %     t_wait = params.t_turn;
    % else 
    %     t_wait = stage.t_ign;
    % end
    t_wait = 0;

    % Retrieve data used multiple times 
    t_burn_tot = stage.t_burn_tot;
    Re = params.Re;

    % Mass estimation
    if t <= t_wait
        m = stage.m0;
    elseif t > t_wait && t <= t_burn_tot + t_wait
        m = stage.m0 - stage.m_dot * (t-t_wait);
    else
        m = stage.m0 - stage.m_dot * t_burn_tot;
    end
    
    % Environment data
    g = params.g0/((1+z/Re)^2);                     % [m/s^2]   - Gravity acceleration taking into account altitude
    rho = getDensity(z);                            % [kg/m^3]  - Density at current altitude
    qdyn = 0.5*rho*velsNorm^2;                      % [Pa]      - Dynamic pressure
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

    if z < 50e3
        angle = 90;
    else
        angle = 45;
    end
    corrector = abs(theta - deg2rad(angle));
    if theta < deg2rad(angle)
        delta = deg2rad(3*corrector);
    else
        delta = deg2rad(-3*corrector);
    end

    % Thrust
    if t > t_wait && t <= t_burn_tot + t_wait       % [N]       - Thrust force acting on the rocket
        T = stage.Thrust;
    else
        T = 0;
    end

    % Forces on the rocket in inertial frame
    F_z = -g*m -D*cos(pi/2-gamma) +L*sin(pi/2-gamma) +T*cos(delta)*sin(theta);
    F_x = -D*sin(pi/2-gamma) -L*cos(pi/2-gamma) +T*cos(delta)*cos(theta);

    % Moment on the rocket
    M_t = T*sin(delta)*(stage.length - stage.xcg) +D*sin(alpha)*(stage.xcp - stage.xcg) -L*cos(alpha)*(stage.xcp - stage.xcg);

    % Derivative vector
    dY = zeros(6, 1);

    dY(1) = xDot;
    dY(2) = zDot;
    dY(3) = F_x/m;
    dY(4) = F_z/m;
    dY(5) = thetaDot;
    dY(6) = M_t/stage.I;

    % Post processing:
    % Forces in body frame
    Fxz_body = dcm'*[F_x F_z]';

    % Prepare output struct for ode recall
    if nargout > 1
        parout.qdyn = qdyn;
        parout.acc = reshape(Fxz_body, [1 2])/m;
        parout.alpha = alpha;
        parout.moment = M_t;
        if exist("t_turn", 'var') && ~isnan(t_turn)
            parout.t_turn = t_turn;
        end
    end
end

