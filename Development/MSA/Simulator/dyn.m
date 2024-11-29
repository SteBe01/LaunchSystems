function [dY, parout] = dyn(t,y, stage, params, current_stage)

    % Y(1) = x (horizontal position, inertial)
    % Y(2) = z (vertical position, inertial)
    % Y(3) = xDot (horizontal velocity, inertial)
    % Y(4) = zDot (vertical velocity, inertial)
    % Y(5) = theta (pitch angle, positive counterclockwise)
    % Y(6) = thetaDot (pitch angular velocity, positive counterclockwise)
    % Y(7) = mass_prop_left

    % Retrieve data from ode
    % x = y(1);
    z = y(2);
    xDot = y(3);
    zDot = y(4);
    theta = y(5);
    thetaDot = y(6);
    m_prop_left = y(7);

    m_prop_left = m_prop_left*(m_prop_left >= 0);

    % Persistent states
    persistent MECO SEPARATION SECO
    if isempty(MECO)
        MECO = false;
    end
    if isempty(SEPARATION)
        SEPARATION = false;
    end
    if isempty(SECO)
        SECO = false;
    end

    if current_stage == 2 && ~SEPARATION && nargout == 1
        fprintf("[%3.1f km] - Stage separation\n", z*1e-3)
        SEPARATION = true;
    end

    velsNorm = norm([xDot zDot]);
    rot_angle = theta - pi/2;
    dcm = [cos(rot_angle) -sin(rot_angle); sin(rot_angle) cos(rot_angle)];
    gamma = atan2(zDot, xDot);
    alpha = -theta+gamma;

    if current_stage == 1
        t_wait = params.t_turn;
    else 
        t_wait = stage.t_ign;
    end

    % Retrieve data used multiple times 
    Re = params.Re;

    % Environment data
    g = params.g0/((1+z/Re)^2);                     % [m/s^2]   - Gravity acceleration taking into account altitude
    rho = getDensity(z);                            % [kg/m^3]  - Density at current altitude
    qdyn = 0.5*rho*velsNorm^2;                      % [Pa]      - Dynamic pressure
    S = pi*(stage.d^2/4);                           % [m^2]     - Rocket surface area
    P = getPressure(z);
    
    % Aerodynamic forces
    D = qdyn*S*stage.Cd;                            % [N]       - Drag force acting on the rocket
    L = qdyn*S*stage.Cl;                            % [N]       - Lift force acting on the rocket
    
    % PID controller
    if current_stage == 1
        if y(2) < 50e3
            angle = 60;
        else
            angle = 45;
        end
    elseif current_stage == 2
        if y(2) < 300e3
            angle = 45;
        else
            angle = 0;
        end
    end
    err = theta - deg2rad(angle);
    delta = -stage.k1*err - stage.k2*thetaDot - stage.k3*alpha;
    if abs(delta) > stage.deltaMax 
        delta = stage.deltaMax*sign(delta);
    end

    % Thrust & mass estimation
    throttling = 1;
    Thrust = interp1(stage.throttling, stage.Thrust, throttling, 'linear', 'extrap');
    m_dot = interp1(stage.throttling, stage.m_dot, throttling, 'linear', 'extrap');
    Pe = interp1(stage.throttling, stage.Pe, throttling, 'linear', 'extrap');

    if t > t_wait && m_prop_left > 0
        T = (Thrust + stage.A_eng*(Pe-P))*stage.N_mot;
    else
        T = 0;
        m_dot = 0;
    end
    m = stage.m0 - stage.m_prop + m_prop_left;

    % Callouts
    if m_prop_left == 0
        if current_stage == 1 && ~MECO && nargout == 1
            fprintf("[%3.1f km] - MECO\n", z*1e-3)
            MECO = true;
        elseif current_stage == 2 && ~SECO && nargout == 1
            fprintf("[%3.1f km] - SECO\n", z*1e-3)
            SECO = true;
        end
    end

    % Compute geom properties
    xcg = interp1(stage.STR_mat(:,1), stage.STR_mat(:,2), m_prop_left, "linear", "extrap");
    I_mat = interp1(stage.STR_mat(:,1), stage.STR_mat(:, 3:5), m_prop_left, "linear", "extrap");
    I = I_mat(2);

    % Forces on the rocket in inertial frame
    F_z = -g*m -D*cos(pi/2-gamma) +L*sin(pi/2-gamma) +T*cos(delta)*sin(theta);
    F_x = -D*sin(pi/2-gamma) -L*cos(pi/2-gamma) +T*cos(delta)*cos(theta);

    % Moment on the rocket
    M_t = T*sin(delta)*(stage.length - xcg) +D*sin(alpha)*(stage.xcp - xcg) + L*cos(alpha)*(stage.xcp - xcg);

    % Derivative vector
    dY = zeros(7, 1);

    dY(1) = xDot;
    dY(2) = zDot;
    dY(3) = F_x/m;
    dY(4) = F_z/m;
    dY(5) = thetaDot;
    dY(6) = M_t/I;
    dY(7) = -m_dot*stage.N_mot;

    % Post processing:
    % Forces in body frame
    Fxz_body = dcm'*[F_x F_z]';

    % Prepare output struct for ode recall
    if nargout > 1
        parout.qdyn = qdyn;
        parout.acc = reshape(Fxz_body, [1 2])/m;
        parout.alpha = alpha;
        parout.moment = M_t;
        parout.rho = rho;
        parout.velssqq = velsNorm^2;
        parout.m = m;
        parout.dv_grav = -g*abs(sin(theta));
        parout.dv_drag = -0.5*S*stage.Cd/stage.m0*rho*velsNorm^2*stage.m0/m;
        parout.delta = delta;
    end
end

