function [dY, parout] = dyn(t,y, stage, params, current_stage)

    % Y(1) = x (horizontal position, inertial)
    % Y(2) = z (vertical position, inertial)
    % Y(3) = xDot (horizontal velocity, inertial)
    % Y(4) = zDot (vertical velocity, inertial)
    % Y(5) = theta (pitch angle, positive counterclockwise)
    % Y(6) = thetaDot (pitch angular velocity, positive counterclockwise)
    % Y(7) = mass_prop_left

    % Retrieve data from ode
    x = y(1);
    z = y(2);
    xDot = y(3);
    zDot = y(4);
    theta = y(5);
    thetaDot = y(6);
    m_prop_left = y(7);

    h = norm([x z]);

    % Retrieve data used multiple times 
    Re = params.Re;

    % Angle of radius vector
    beta = atan2(z, x);

    m_prop_left = m_prop_left*(m_prop_left >= 0);

    % Persistent states
    persistent MECO SEPARATION SECO APOGEE
    if isempty(MECO)
        MECO = false;
    end
    if isempty(SEPARATION)
        SEPARATION = false;
    end
    if isempty(SECO)
        SECO = false;
    end
    if isempty(APOGEE)
        APOGEE = false;
    end

    if current_stage == 2 && ~SEPARATION && nargout == 1
        if params.dispStat
            fprintf("[%3.1f km] - Stage separation\n", (z - Re)*1e-3)
        end
        SEPARATION = true;
    end
    if current_stage == 2 && y(4) < 0 && ~APOGEE && nargout == 1
        if params.dispStat
            fprintf("[%3.1f km, vx = %4.3f km/s] - Apogee\n", (h - Re)*1e-3, y(3)*1e-3)
        end
        APOGEE = true;
    end

    velsNorm = norm([xDot zDot]);
    rot_angle = theta;
    dcm = [cos(rot_angle) -sin(rot_angle); sin(rot_angle) cos(rot_angle)];
    gamma = atan2(zDot, xDot);
    alpha = theta-gamma;

    if current_stage == 1
        t_wait = params.t_turn;
    else 
        t_wait = stage.t_ign;
    end

    % Environment data
    [~, a, P, rho] = computeAtmosphericData(h - Re);
    g = 398600*1e9/h^2;
    if abs(g) > 1e2
        warning("Module of gravity has gone out of safe values. g: " + num2str(g));
    end

    qdyn = 0.5*rho*velsNorm^2;                      % [Pa]      - Dynamic pressure
    S = pi*(stage.d^2/4);                           % [m^2]     - Rocket surface area
    M = velsNorm/a;

    % Compute all necessary interpolations
    
    % Geometric parameters
    xcg = interp1(stage.STR_mat(:,1), stage.STR_mat(:,2), m_prop_left, "linear", "extrap");
    I_mat = interp1(stage.STR_mat(:,1), stage.STR_mat(:, 3:5), m_prop_left, "linear", "extrap");
    I = I_mat(2);

    % Engine parameters 
    throttling = 1;
    Thrust = interp1(stage.throttling, stage.Thrust, throttling, 'linear', 'extrap');
    m_dot = interp1(stage.throttling, stage.m_dot, throttling, 'linear', 'extrap');
    Pe = interp1(stage.throttling, stage.Pe, throttling, 'linear', 'extrap');

    % Aerodynamic coefficients
    if current_stage == 1
        interpValues = params.coeffs({1:4, M, rad2deg(alpha), (h - Re)});
        Cd = interpValues(1)*params.CD_mult;
        Cl = interpValues(2)*params.CL_mult;
        xcp = interpValues(3);
    else
        Cd = 0;
        Cl = 0;
        xcp = 4;
    end

    % Aerodynamic forces
    D = qdyn*S*Cd;                            % [N]       - Drag force acting on the rocket
    L = qdyn*S*Cl;                            % [N]       - Lift force acting on the rocket

    % Thrust & mass estimation
    if t > t_wait && m_prop_left > stage.m_prop_final
        T = (Thrust + stage.A_eng*(Pe-P))*stage.N_mot;
    else
        T = 0;
        m_dot = 0;
    end
    m = stage.m0 - stage.m_prop + m_prop_left;

    % Callouts
    if m_prop_left <= stage.m_prop_final
        if current_stage == 1 && ~MECO && nargout == 1
            if params.dispStat
                fprintf("[%3.1f km] - MECO\n", (z - Re)*1e-3)
            end
            MECO = true;
        elseif current_stage == 2 && ~SECO && nargout == 1
            if params.dispStat
                fprintf("[%3.1f km] - SECO\n", (z - Re)*1e-3)
            end
            SECO = true;
        end
    end
    
    % PID controller
    if T ~= 0
        angle = getPitch(params.pitch, (h - Re));
        err = theta - angle;
        delta = -stage.k1*err - stage.k2*thetaDot - stage.k3*alpha;
        delta = -delta;
        if abs(delta) > stage.deltaMax
            delta = stage.deltaMax*sign(delta);
        end
    else
        delta = 0;
    end

    % Forces on the rocket in inertial frame
    F_x = L*sin(alpha)*sign(alpha) - D*cos(alpha) + T*cos(delta);
    F_z = L*cos(alpha)*sign(alpha) + D*sin(alpha) + T*sin(alpha);
    F_body = [F_x F_z]';
    F_in = dcm*F_body;
    F_x = F_in(1) - m*g*cos(beta);
    F_z = F_in(2) - m*g*sin(beta);

    F_L_in = dcm*[L*sin(alpha)*sign(alpha) L*cos(alpha)*sign(alpha)]';
    F_D_in = dcm*[-D*cos(alpha) D*sin(alpha)]';

    % Moment on the rocket
    M_t = - L*cos(alpha)*(xcp - xcg) - D*sin(alpha)*(xcp - xcg) - T*sin(delta)*(stage.length - xcg);

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
        parout.dv_drag = -0.5*S*Cd/stage.m0*rho*velsNorm^2*stage.m0/m;
        parout.delta = delta;
        parout.coeffs = [Cd Cl xcp];

        parout.dcm = dcm;
        parout.F_in = F_in';
        parout.F_L_in = F_L_in;
        parout.F_D_in = F_D_in;
    end
end

