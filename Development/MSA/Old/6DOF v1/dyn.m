function [dY, parout] = dyn(t,Y, stage, params)

    % Y(1) = x (horizontal position, inertial)
    % Y(2) = y (horizontal position, inertial)
    % Y(3) = z (vertical position, inertial)
    % Y(4) = xDot (horizontal velocity, inertial)
    % Y(5) = yDot (horizontal velocity, inertial)
    % Y(6) = zDot (vertical velocity, inertial)
    % Y(7) = x_theta (pitch angle, positive counterclockwise)
    % Y(8) = y_theta (pitch angle, positive counterclockwise)
    % Y(9) = z_theta (pitch angle, positive counterclockwise)
    % Y(10) = x_thetaDot (pitch angular velocity, positive counterclockwise)
    % Y(11) = y_thetaDot (pitch angular velocity, positive counterclockwise)
    % Y(12) = z_thetaDot (pitch angular velocity, positive counterclockwise)

    % Retrieve data from ode
    x = Y(1);
    y = Y(2);
    z = Y(3);
    xDot = Y(4);
    yDot = Y(5);
    zDot = Y(6);
    x_theta = Y(7);
    y_theta = Y(8);
    z_theta = Y(9);
    x_thetaDot = Y(10);
    y_thetaDot = Y(11);
    z_thetaDot = Y(12);

    velsNorm = norm([xDot yDot zDot]);
    
    R1 = [1 0 0
        0 cos(x_theta) -sin(x_theta)
        0 sin(x_theta) cos(x_theta)];
    R2 = [cos(y_theta) 0 sin(y_theta)
        0 1 0
        -sin(y_theta) 0 cos(y_theta)];
    R3 = [cos(z_theta) -sin(z_theta) 0
        sin(z_theta) cos(z_theta) 0
        0 0 1];
    R_in_body = R3*R2*R1;
    % R_in_body = R_in_body';
    R_body_in = R_in_body';

    gamma1_in = rad2deg(atan2(xDot,yDot));
    gamma2_in = rad2deg(atan2(-zDot,sqrt(yDot^2+xDot^2)));

    V_in = [xDot yDot zDot]';
    V_body = R_in_body*V_in;

    gamma1_body = atan2(V_body(1),V_body(2));
    gamma2_body = -atan2(V_body(3),sqrt(V_body(2)^2+V_body(1)^2));

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
    g = params.g0/((1-z/Re)^2);                     % [m/s^2]   - Gravity acceleration taking into account altitude
    rho = getDensity(-z);                            % [kg/m^3]  - Density at current altitude
    qdyn = 0.5*rho*velsNorm^2;                      % [Pa]      - Dynamic pressure
    S = pi*(stage.d^2/4);                           % [m^2]     - Rocket surface area
    
    % Aerodynamic forces
    D = qdyn*S*stage.Cd;                            % [N]       - Drag force acting on the rocket
    L = qdyn*S*stage.Cl;                            % [N]       - Lift force acting on the rocket
    
    % Thrust
    if t <= 0
        T = stage.Thrust;
    else
        T = 0;
    end

    delta2 = deg2rad(0);
    delta1 = deg2rad(0);

    b1 = 10;
    b2 = 6;

    g_in = g*[0 0 1]';
    g_body = R_in_body*g_in;

    % Forces on the rocket in body frame
    F_x_1 = T*cos(delta2)*cos(delta1);
    F_y_1 = -T*cos(delta2)*sin(delta1);
    F_z_1 = T*sin(delta2);
    % Moment on the rocket
    % M_x_1 = (stage.J - stage.K)*z_thetaDot*y_thetaDot;
    % M_y_1 = F_z_1*b1 + (stage.K - stage.I)*x_thetaDot*z_thetaDot;
    % M_z_1 = -F_y_1*b1 + (stage.I - stage.J)*x_thetaDot*y_thetaDot;
    M_x_1 = 0;
    M_y_1 = F_z_1*b1;
    M_z_1 = -F_y_1*b1;
    % Forces on the rocket in body frame
    F_x_2 = -D*cos(gamma2_body)*sin(gamma1_body) -L*sin(gamma2_body)*cos(gamma1_body);
    F_y_2 = -D*cos(gamma2_body)*cos(gamma1_body) -L*sin(gamma2_body)*sin(gamma1_body);
    F_z_2 = D*sin(gamma2_body) -L*cos(gamma2_body);
    % Moment on the rocket
    M_x_2 = 0;
    M_y_2 = F_z_2*b2;
    M_z_2 = -F_y_2*b2;

    F_x = F_x_1 + F_x_2 + g_body(1)*m;
    F_y = F_y_1 + F_y_2 + g_body(2)*m;
    F_z = F_z_1 + F_z_2 + g_body(3)*m;
    M_x = M_x_1 + M_x_2;
    M_y = M_y_1 + M_y_2;
    M_z = M_z_1 + M_z_2;

    F_in = R_body_in*[F_x F_y F_z]';
    F_x = F_in(1);
    F_y = F_in(2);
    F_z = F_in(3);
    M_in = R_body_in*[M_x M_y M_z]';
    M_x = M_in(1);
    M_y = M_in(2);
    M_z = M_in(3);

    % Derivative vector
    dY = zeros(12, 1);

    dY(1) = xDot;
    dY(2) = yDot;
    dY(3) = zDot;
    dY(4) = F_x/m;
    dY(5) = F_y/m;
    dY(6) = F_z/m;
    dY(7) = x_thetaDot;
    dY(8) = y_thetaDot;
    dY(9) = z_thetaDot;
    dY(10) = M_x/stage.I;
    dY(11) = M_y/stage.J;
    dY(12) = M_z/stage.K;

    % Post processing:
    % Forces in body frame
    % Fxz_body = dcm'*[F_x F_z]';

    % Prepare output struct for ode recall
    if nargout > 1
        parout.dcm = R_in_body;
    end
end

