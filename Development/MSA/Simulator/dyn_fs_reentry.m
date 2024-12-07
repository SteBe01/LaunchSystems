function [dY, parout] = dyn_fs_reentry(t,y, stages, params, varargin)
% VALID ONLY FOR THE FIRST STAGE!

    % Frame of reference: inertial with origin in Earth center
    % Y(1) = x (horizontal position, inertial)                              [m]
    % Y(2) = z (vertical position, inertial)                                [m]
    % Y(3) = xDot (horizontal velocity, inertial)                           [m/s]
    % Y(4) = zDot (vertical velocity, inertial)                             [m/s]
    % Y(5) = theta (pitch angle, positive counterclockwise)                 [rad]
    % Y(6) = thetaDot (pitch angular velocity, positive counterclockwise)   [rad]
    % Y(7) = mass_prop_left (propellant mass left in the tanks)             [kg]
    % Y(8) = temperature                                                    [k]
    stage = stages.stg1;

    % Retrieve data from ode
    x = y(1);
    z = y(2);
    xDot = y(3);
    zDot = y(4);
    theta = y(5);
    thetaDot = y(6);
    m_prop_left = y(7);
    Temp = y(8);
    % If propellant mass is negative, set it to zero
    m_prop_left = m_prop_left*(m_prop_left >= 0); 

    % Pre-define derivative vector
    dY = zeros(8,1);

    % Radius [m]
    h = norm([x z]);

    % Earth radius [m] 
    Re = params.Re;

    % Angle of radius vector [rad]
    beta = atan2(z, x);

    % Angle used for PID computation [rad]
    xi = theta + pi/2 - beta;

    % Velocities in the rotating frame
    ang = pi/2-beta;
    vec_rotated = [cos(ang) -sin(ang); sin(ang) cos(ang)]*[xDot; zDot];

    % Gravitational attraction [m/s^2]
    g = 398600*1e9/h^2;

    % Velocities norm [m/s]
    velsNorm = norm([xDot zDot]);

    % dcm from inertial to body
    rot_angle = theta;
    dcm = [cos(rot_angle) -sin(rot_angle); sin(rot_angle) cos(rot_angle)];
    
    % Flight path angle [rad]
    gamma = atan2(zDot, xDot);
    % Angle wrt relative velocity [rad]
    alpha = theta-gamma;

    % Environment data
    [Tatm, a, P, rho] = computeAtmosphericData(h - Re);
    % Mach number [-]
    M = velsNorm/a;

    % Dynamic pressure [Pa]
    qdyn = 0.5*rho*velsNorm^2;
    % Rocket surface area [m^2]
    S = pi*(stage.d^2/4);

    % Compute all necessary interpolations
    
    % Geometric parameters wrt propellant mass left:
    xcg = 6;
    I = 1e3;

    % Engine parameters wrt throttling percentage
    throttling = stage.prcg_throt;
    Thrust = interp1(stage.throttling, stage.Thrust, throttling, 'linear', 'extrap');
    m_dot_interp = interp1(stage.throttling, stage.m_dot, throttling, 'linear', 'extrap');
    Pe = interp1(stage.throttling, stage.Pe, throttling, 'linear', 'extrap');

    % Aerodynamic coefficients
    interpValues = params.coeffs({1:4, M, rad2deg(alpha), (h - Re)});
    Cd = interpValues(1)*params.CD_mult;
    Cl = interpValues(2)*params.CL_mult;
    Cd = 1;
    Cl = 0;
    % xcp = interpValues(3);
    xcp = 4;

    margin = xcp - xcg;
    if abs(alpha) < deg2rad(90)
        margin = -margin;
    end

    % Aerodynamic forces
    D = qdyn*S*Cd;                            % [N]       - Drag force acting on the rocket
    L = qdyn*S*Cl;                            % [N]       - Lift force acting on the rocket
    
    h1 = 60e3 + Re;
    h2 = 9e3 + Re;
    h3 = 1e3 + Re;
    Cd1 = 0.6;
    Cd2 = 0.8;
    Cd3 = 0.8;
    S1 = 20;
    S2 = 44.8;
    S3 = 498.8;

    if h < h3
        Cd_par = Cd3;
        S_par = S3;
    elseif h < h2
        Cd_par = Cd2;
        S_par = S2;
    elseif h < h1
        Cd_par = Cd1;
        S_par = S1;
    else
        Cd_par = 0;
        S_par = 0;
    end

    D_par = qdyn*S_par*Cd_par;

    % Thrust & mass estimation
    T = 0;
    m_dot = 0;

    % m = stage.m0 - stage.m_prop + m_prop_left;
    m = (stage.m0-stages.stg2.m0) - stage.m_prop + m_prop_left;
    
    % PID controller
    delta = 0;

    % Forces on the rocket in inertial frame
    F_x = L*sin(alpha)*sign(alpha) - D*cos(alpha) - D_par*cos(alpha) + T*cos(delta);
    F_z = L*cos(alpha)*sign(alpha) + D*sin(alpha) + D_par*sin(alpha) + T*sin(delta);
    F_body = [F_x F_z]';
    F_in = dcm*F_body;
    F_x = F_in(1) - m*g*cos(beta);
    F_z = F_in(2) - m*g*sin(beta);

    F_L_in = dcm*[L*sin(alpha)*sign(alpha) L*cos(alpha)*sign(alpha)]';
    F_D_in = dcm*[-D*cos(alpha) D*sin(alpha)]';

    % Moment on the rocket
    M_t = - L*cos(alpha)*margin - D*sin(alpha)*margin + D_par*sin(alpha)*(stage.length - xcg) - T*sin(delta)*(stage.length - xcg);

    % Heat flux
    Length = stages.stg1.length - stages.stg2.length;
    sigma = 5.670400* 10^-8;                       % Stefan-Boltzmann
    eps = 0.9; % Emissivity of the surface
    cp =  900;  % Specific heat capacity (J/kg/K)
    Qdot_aer = 0.5*rho*velsNorm^3*S;
    Qdot_space = eps*sigma*(Temp^4-Tatm^4)*(2*S + pi*stage.d*Length);
    Qdot_tot = Qdot_aer-Qdot_space;
    b = 5.5164* 10^-5; 
    q_dot_stagn = b * sqrt(rho)/sqrt(stage.d)*velsNorm^3.15;

    % Derivative vector

    dY(1) = xDot;
    dY(2) = zDot;
    dY(3) = F_x/m;
    dY(4) = F_z/m;
    dY(5) = thetaDot;
    dY(6) = M_t/I;
    dY(7) = -m_dot*stage.N_mot;
    dY(8) = Qdot_tot/(m*cp);
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

        dcm_xi = [cos(xi) -sin(xi); sin(xi) cos(xi)];
        parout.dcm = dcm_xi;
        parout.F_in = F_in';
        parout.F_L_in = F_L_in;
        parout.F_D_in = F_D_in;
        parout.Thrust = T;
        % parout.Q = Qdot_aer/S + Qdot_space/(2*S + pi*stage.d*Length);
        parout.Q = Qdot_tot;
    end
end

