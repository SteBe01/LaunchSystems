function [dY, parout] = dyn(t,y, stage, params, current_stage, varargin)

    % Frame of reference: inertial with origin in Earth center
    % Y(1) = x (horizontal position, inertial)                              [m]
    % Y(2) = z (vertical position, inertial)                                [m]
    % Y(3) = xDot (horizontal velocity, inertial)                           [m/s]
    % Y(4) = zDot (vertical velocity, inertial)                             [m/s]
    % Y(5) = theta (pitch angle, positive counterclockwise)                 [rad]
    % Y(6) = thetaDot (pitch angular velocity, positive counterclockwise)   [rad]
    % Y(7) = mass_prop_left (propellant mass left in the tanks)             [kg]

    % Retrieve data from ode
    x = y(1);
    z = y(2);
    xDot = y(3);
    zDot = y(4);
    theta = y(5);
    thetaDot = y(6);
    m_prop_left = y(7);
    % If propellant mass is negative, set it to zero
    m_prop_left = m_prop_left*(m_prop_left >= 0); 

    % Pre-define derivative vector
    dY = zeros(7,1);

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

    % Persistent states
    persistent MECO SEPARATION SECO APOGEE FRANCO END_BURN
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

    if isempty(FRANCO)
        FRANCO = false;
    end
    if isempty(END_BURN)
        END_BURN = false;
    end

    % Callouts
    if current_stage == 2 && ~SEPARATION && nargout == 1
        if params.dispStat
            fprintf("[%3.1f km] - Stage separation\n", (h - Re)*1e-3)
        end
        SEPARATION = true;
    end
    if current_stage == 2 && vec_rotated(2) < 0 && ~APOGEE && nargout == 1
        if params.dispStat
            fprintf("[%3.1f km, vx = %4.3f km/s] - Apogee\n", (h - Re)*1e-3, vec_rotated(1)*1e-3)
        end
        APOGEE = true;
    end

    % Gravitational attraction [m/s^2]
    g = 398600*1e9/h^2;

    % Fast intertial orbit propagation
    % The last part of inertial flight can be computed without requiring
    % all the other computations, in order to save computational time
    if nargout == 1 && END_BURN
        dY(1) = xDot;
        dY(2) = zDot;
        dY(3) = -g*cos(beta);
        dY(4) = -g*sin(beta);
        dY(5) = thetaDot;
        dY(6) = 0;
        dY(7) = 0;

        return
    end

    % Velocities norm [m/s]
    velsNorm = norm([xDot zDot]);

    % dcm from inertial to body
    rot_angle = theta;
    dcm = [cos(rot_angle) -sin(rot_angle); sin(rot_angle) cos(rot_angle)];
    
    % Flight path angle [rad]
    gamma = atan2(zDot, xDot);
    % Angle wrt relative velocity [rad]
    alpha = theta-gamma;

    % Wait time before ignition:
    % - For the first stage it's the time between detach from airplane and ignition
    % - For the second stage it's the time between stage separation and ignition
    if current_stage == 1
        t_wait = params.t_turn;
    else 
        t_wait = stage.t_ign;
    end

    % Environment data
    % a: Speed of sound [m/s]
    % P: Pressure [Pa]
    % rho: Density [kg/m^3]
    [~, a, P, rho] = computeAtmosphericData(h - Re);
    % Mach number [-]
    M = velsNorm/a;

    % Dynamic pressure [Pa]
    qdyn = 0.5*rho*velsNorm^2;
    % Rocket surface area [m^2]
    S = pi*(stage.d^2/4);

    % Compute all necessary interpolations
    
    % Geometric parameters wrt propellant mass left:
    % xcg: position of center of mass (starting from the nosecose) [m]
    % I_mat: Matrix of inertia moments [kg*m^2]
    xcg = interp1(stage.STR_mat(:,1), stage.STR_mat(:,2), m_prop_left, "linear", "extrap");
    I_mat = interp1(stage.STR_mat(:,1), stage.STR_mat(:, 3:5), m_prop_left, "linear", "extrap");
    I = I_mat(2);

    % Engine parameters wrt throttling percentage
    % Thrust: Thrust without considering expansion [N]
    % m_dot_interp: Mass flow rate [kg/s]
    % Pe: Nozzle exit pressure [Pa]
    throttling = stage.prcg_throt;
    if h-Re > params.h_reign && abs(xi) > params.xi_err 
        throttling = 0.1;
    end
    Thrust = interp1(stage.throttling, stage.Thrust, throttling, 'linear', 'extrap');
    m_dot_interp = interp1(stage.throttling, stage.m_dot, throttling, 'linear', 'extrap');
    Pe = interp1(stage.throttling, stage.Pe, throttling, 'linear', 'extrap');

    % Aerodynamic coefficients
    % Cd: Drag coefficient [-]
    % Cl: Lift coefficient [-]
    % xcp: Position of center of pressure (starting from nosecone) [m]
    if current_stage == 1 && t > t_wait
        interpValues = params.coeffs({1:4, M, rad2deg(alpha), (h - Re)});
        Cd = interpValues(1)*params.CD_mult;
        Cl = interpValues(2)*params.CL_mult;
        % xcp = interpValues(3);
        xcp = 6;
    else
        Cd = 0;
        Cl = 0;
        xcp = 4;
    end

    % Aerodynamic forces
    D = qdyn*S*Cd;                            % [N]       - Drag force acting on the rocket
    L = qdyn*S*Cl;                            % [N]       - Lift force acting on the rocket

    % Thrust & mass estimation
    if t > t_wait && m_prop_left > stage.m_prop_final && not(current_stage == 2 && h-Re > params.h_shutoff)
        T = (Thrust + stage.A_eng*(Pe-P))*stage.N_mot;
        m_dot = m_dot_interp;
    else
        T = 0;
        m_dot = 0;
    end

    v_orb = sqrt(398600/(Re*1e-3+400))*1e3;
    v_thr = 1;
    if nargout > 1
        v_thr = 3;
    end
    if h-Re > params.h_reign && abs(vec_rotated(1) - v_orb) > v_thr && ~END_BURN
        FRANCO = true;
        T = (Thrust + stage.A_eng*(Pe-P))*stage.N_mot;
        m_dot = m_dot_interp;
    else
        if FRANCO == true && h-Re > params.h_reign && ~END_BURN
            END_BURN = true;
        end
    end

    if current_stage == 2 && (m_prop_left <= 0.05*stage.m_prop)
        T = 0;
        m_dot = 0;
    end
    if END_BURN
        T = 0;
        m_dot = 0;
    end

    m = stage.m0 - stage.m_prop + m_prop_left;

    % Callouts
    if m_prop_left <= stage.m_prop_final
        if current_stage == 1 && ~MECO && nargout == 1
            if params.dispStat
                fprintf("[%3.1f km] - MECO\n", (h - Re)*1e-3)
            end
            MECO = true;
        elseif current_stage == 2 && ~SECO && nargout == 1
            if params.dispStat
                fprintf("[%3.1f km] - SECO\n", (h - Re)*1e-3)
            end
            SECO = true;
        end
    end
    
    % PID controller
    if T ~= 0
        if ~isempty(varargin)
            angle = varargin{1};
        else
            angle = getPitch(params.pitch, (h - Re));
        end
        if h-Re > params.h_reign
            angle = 0;
        end
        err = xi - angle;
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
    F_z = L*cos(alpha)*sign(alpha) + D*sin(alpha) + T*sin(delta);
    F_body = [F_x F_z]';
    F_in = dcm*F_body;
    F_x = F_in(1) - m*g*cos(beta);
    F_z = F_in(2) - m*g*sin(beta);

    F_L_in = dcm*[L*sin(alpha)*sign(alpha) L*cos(alpha)*sign(alpha)]';
    F_D_in = dcm*[-D*cos(alpha) D*sin(alpha)]';

    % Moment on the rocket
    M_t = - L*cos(alpha)*(xcp - xcg) - D*sin(alpha)*(xcp - xcg) - T*sin(delta)*(stage.length - xcg);

    % Derivative vector

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

        dcm_xi = [cos(xi) -sin(xi); sin(xi) cos(xi)];
        parout.dcm = dcm_xi;
        parout.F_in = F_in';
        parout.F_L_in = F_L_in;
        parout.F_D_in = F_D_in;
        parout.Thrust = T;
    end
end

