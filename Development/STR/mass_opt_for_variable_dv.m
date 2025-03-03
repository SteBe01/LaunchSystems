

%% the iterative method verifies mass estimation and assesses optimal staging

clc;
clear;
close all;

%data from other departments:
Is = [324.5; 343]; %[s] stages Is
M.pay = 250; %[kg] nominal payload mass
M.pay_max = 400; %[kg] maximum payload mass
M.adapter = 0.0755 * M.pay + 50; %[kg] estimated mass from Edberg-Costa
M.pay_effective = M.pay + M.adapter; %[kg] Edberg-Costa includes adapter mass in the payload
% M.pay_effective = M.pay; %[kg] Edberg-Costa includes adapter mass in the payload
OF = 2.58; %[-] Ox/Fu ratio for LOX-RP1
loads.nx = 6;%8;%6;  %longitudinal acceleration load factor [-]
loads.nz = 0.8; %1.2;%1; %transversal acceleration load factor [-]
loads.K = 1.50; %loads resistance safety factor [-]
FoS1 = 1.50; %FoS of first stage [-]
FoS2 = 1.50; %FoS of first stage [-]
S_tail = 2; %[m^2]
S_ft2 = S_tail * 10.7639; %[ft^2]
m_tail = (5*0.453592) * S_ft2^1.09; %[kg] mass of the tail of the rocket

a = 0.5; %0.71; %balancing factor for iterative method
n = 1;
m = 100;
dv_it = linspace(9.725, 9.725, n);
diam1_it = linspace(1, 1.6, m);
Nmax = 200; %maximum number of iterations for the while loop

% M_it = zeros(n*m,1);
% h_it = zeros(n*m,1);
% loads_it = zeros(n*m,1);

for j = 1:n
    for k = 1:m
        skipcycle = false;
        dv = dv_it(j); %[km/s] required dv
        diam1 = diam1_it(k); %[m] external diameter of first stage
        diam2 = 0.95 * diam1; %[m] external diameter of second stage and fairing
        AR = sqrt(2);%sqrt(3); %aspect ratio of oblate domes [-]
        maxQ = 50000; %0.5 * rho_air * v^2; %[Pa] maximum dynamic pressure
        Ca1 = 1.3; %drag coefficient of first stage
        Ca2 = 1.3; %drag coefficient of second stage
        Caf = 1.3; %drag coefficient of the fairing
        S1 = pi * diam1^2 / 4; %first stage cross section [m^2]
        S2 = pi * diam2^2 / 4; %second stage cross section [m^2]
        loads.F_drag = zeros(3, 1); %aerodynamic forces acting on the three part of the launcher (fairing, stg2, stg1) [N]
        TW1 = 1.3; %[-] T/W ratio of first stage
        TW2 = 0.8; %[-] T/W ratio of second stage
        T1 = 27.4 * 1e3; %[N] 1 electron - rutherford motor thrust (first stage)
        T2 = 31 * 1e3; %[N] 1 electron - rutherford motor thrust (second stage)
        
        %stage 1
        M1.OF = OF;%[-] Ox/Fu ratio
        %M1.motor = 450; %[kg] only motor, pumps and batteries (electron - rutherford motor) %pump-fed
        M1.m_dot = 7.66; %[kg/s] mass propellant flow rate of each motor
        M1.rhorp1 = 820;  %[kg/m^3] density of rp1
        M1.rholox = 1140; %[kg/m^3] density of lox
        M1.avionics = 75 * 0.2; %[kg] from Edberg-Costa
        M1.other = 0; %185.7; %3 * 64.2; %256.8; %218; %[kg] 
        M1.wing = m_tail; %[kg]
        M1.stg = 1; %[#] stage ID
        h1.motor = 0.75; %[m] height of the motor
        h1.h0 = 0; %[m] starting height
        mat1 = 10; % 1 for Ti, 2 for Al 2XXX, 3 for Steel, 4 for Carbon Fiber Toray M46J, 5 for Al 7075 T6, 6 for Al 2090, 7 for CF Hexcel® HexTow® IM7, 8 for Al 6061 T6, 9 for 300M Steel alloy, 10 for Al 2219, 11 for Carbon/Epoxy (MatWeb) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% in future versions can be optimized the material selection in function
        press1 = 2; % 0 for unpressurized, 1 for pressure-fed, 2 for pump-fed, 3 for blowdown
        
        %stage 2
        M2.OF = OF;%[-] Ox/Fu ratio
        %M2.motor = 45; %[kg] only motor, pumps and batteries (electron - rutherford motor) %pump-fed
        M2.m_dot = 7.66; %[kg/s] mass propellant flow rate of each motor
        M2.rhorp1 = 820;  %[kg/m^3] density of rp1
        M2.rholox = 1140; %[kg/m^3] density of lox
        M2.avionics = 75 * 0.8; %[kg] from Edberg-Costa
        M2.other = 0; %[kg] 
        M2.wing = 0; %[m^2]
        M2.stg = 2; %[#] stage ID
        h2.motor = 0.89; %[m] height of the motor
        h2.h0 = 12; %[m] starting height
        mat2 = 11; % 1 for Ti, 2 for Al 2XXX, 3 for Steel, 4 for Carbon Fiber Toray M46J, 5 for Al 7075 T6, 6 for Al 2090, 7 for CF Hexcel® HexTow® IM7, 8 for Al 6061 T6, 9 for 300M Steel alloy, 10 for Al 2219, 11 for Carbon/Epoxy (MatWeb)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% in future versions can be optimized the material selection in function
        press2 = 2; % 0 for unpressurized, 1 for pressure-fed, 2 for pump-fed, 3 for blowdown
        
        %fairing:
        fairing.mat_id = 4; % 1 for Ti, 2 for Al 2XXX, 3 for Steel, 4 for Carbon Fiber, 5 for Al 7XXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% in future versions can be optimized the material selection in function
        
        %first guesses:
        eps0 = [0.10; 0.1]; %[0.06; 0.2]; %[-] stages structural mass indexes
        fairing.base_diam = diam2; %[m] first guess for fairing base diameter
        M.M0 = 17000; %[kg]
        M.M1 = 1800; %[kg]
        h.tot = 20; %[m]
        h.CG = 9; %[m]
        h1.tot = 14; %[m]
        h2.tot = 5; %[m]
        h1.attach = 15; %[m]
        h2.attach = 19; %[m]
        
        %while loop parameters:
        i = 2;
        err = 1;
        tol = 1e-8;
        eps_real = zeros(2, Nmax-1);
        eps_real(:,1) = eps0;
        
        while i < Nmax && err > tol
        %optimal staging:
        [m_stag, m_init, m_prop] = tandem_opt_staging(Is, eps_real(:, i-1), dv, M.pay_effective, 0);

        if m_init > 1e10
            %M = 0;
            %h = 0;
            skipcycle = true;
            break
        end

        %compute stage parameters:
        M1.motor = ceil(TW1 * M.M0 * 9.81 / T1) * 40; %[kg] total motor mass of first stage
        M2.motor = ceil(TW2 * M.M1 * 9.81 / T2) * 45; %[kg] total motor mass of second stage
        M1.n_mot = M1.motor / 40; %[-] number of motors of first stage
        M2.n_mot = M2.motor / 45; %[-] number of motors of second stage
        M1.Thrust = M1.n_mot * T1; %[N] total thrust of first stage
        M2.Thrust = M2.n_mot * T2; %[N] total thrust of second stage
    
        %recover MER data:
        M.wiring = 1.43 * h.tot; %[kg] from Edberg-Costa
    
        %divide MER data between stages:
        M1.wiring = M.wiring * h1.tot / (h1.tot + h2.tot);%[kg]
        M2.wiring = M.wiring * h2.tot / (h1.tot + h2.tot);%[kg]
        M1.T_struct = 2.55 * M1.Thrust * 1e-4; %[kg]
        M2.T_struct = 2.55 * M2.Thrust * 1e-4; %[kg]
    
        %recover GLOM data:
        M1.prop = m_prop(1);
        M2.prop = m_prop(2);
    
        %compute bending loads:
        % loads.M_exp = loads.nz * M.M0 * 9.81 * h.tot / 5;
        loads.M_max = loads.nz * M.M0 * 9.81 * h.CG * ( 1 - h.CG / h.tot );
    
        % %fairing:
        % Sf = pi * diam2^2 / 4; %fairing cross section [m^2]
        % loads.F_drag(3) = maxQ * Caf * Sf; %compose aerodynamic forces vector [N]
        % loads_f.F_drag = loads.F_drag(3); %aerodynamic force acting on fairing [N]
        % [fairing] = fairing_fun(M.pay_max, M.pay, fairing);
    
        %fairing
        Sf = pi * diam2^2 / 4; %fairing cross section [m^2]
        fairing.base_diam = diam2;
        loads_f = loads;
        loads.F_drag(3) = maxQ * Caf * Sf; %compose aerodynamic forces vector [N]
        loads_f.F_drag = loads.F_drag(3); %aerodynamic force acting on fairing [N]
        [fairing, loads_f] = fairing_fun1(M.pay, fairing, loads_f);
    
        %stage 2:
        M2.R_next = diam2 / 2; %[m]
        M2.fairing = 0; %fairing.m; %[kg] WE ASSUME THE ADAPTER DETATCH FROM THE LAUNCHER TOGETHER WITH THE SECOND STAGE
        loads2 = loads; %recover loads
        loads2.K = FoS2; %Factore of safety of second stage [-]
        loads2.m = M2.avionics + M2.wiring + M.pay_effective + fairing.m;%sustained mass [kg]
        % loads2.h_m = fairing.h_m; %[m] barycenter height of sustained mass %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        if h1.attach > h.CG
            loads2.M_exp = loads.M_max * ( h.tot - h1.attach ) / ( h.tot - h.CG );
        else
            loads2.M_exp = loads.M_max * h1.attach / h.CG;
        end
        if (S2-Sf) > 0
            loads.F_drag(2) = maxQ * Ca2 * (S2-Sf); %compose aerodynamic forces vector [N]
        else
            loads.F_drag(2) = 0; %compose aerodynamic forces vector [N]
        end
        loads2.F_drag = sum( loads.F_drag ); %aerodynamic force [N]
        [M2, h2, th2] = inert_mass_common_dome(M2, h2, diam2, AR, loads2, mat2, press2);
        
        if M2.diam_cyl < diam1
            diam2_old = diam2;
            diam2 = M2.diam_cyl; %[m] correct the diameter to have a cylindrical tank on the second stage
        else
            diam2 = diam1;
        end
    
        % if diam2 > diam1
        %     diam1 = diam2;
        % end

        %stage 1:
        M1.R_next = M2.R_end; %[m]
        M1.fairing = fairing.m; %0; %[kg] WE ASSUME THE FAIRING DETATCH FROM THE LAUNCHER TOGETHER WITH THE FIRST STAGE
        h1.motor_next = h2.motor; %[m]
        loads1 = loads; %recover loads
        loads1.K = FoS1; %Factore of safety of second stage [-]
        loads1.m = M1.avionics + M1.wiring + M2.tot + M.pay_effective + M1.fairing;%sustained mass [kg] M2.tot comprises the adapter
        % loads1.h_m = 
        loads1.M_exp = loads.M_max; 
        if (S1-S2) > 0
            loads.F_drag(1) = maxQ * Ca1 * (S1-S2); %compose aerodynamic forces vector [N]
        else
            loads.F_drag(1) = 0; %compose aerodynamic forces vector [N]
        end
        loads1.F_drag = sum( loads.F_drag ); %aerodynamic force [N]
        [M1, h1, th1] = inert_mass_common_dome(M1, h1, diam1, AR, loads1, mat1, press1);
    
        %recover eps_real:
        eps_real(:,i) = (1-a) * eps_real(:,i-1) + (a) * [M1.eps; M2.eps];
    
        %recover real dimensions of M1, M2, of the fairing and of the adapter:
        M.M0 = M.pay_effective + M1.tot + M2.tot;%[kg] initial mass
        M.M1 = M.M0 - M1.tot; %[kg] mass after first stage separation
        M.eps1 = M1.eps; %[-] first stage eps
        M.eps2 = M2.eps; %[-] second stage eps
        M.tb1 = M1.tb; %[s] first stage burning time
        M.tb2 = M2.tb; %[s] second stage burning time
        M.nx1 = M1.n_final; %[-] ending longitudinal load factor
        M.nx2 = M2.n_final; %[-] ending longitudinal load factor
        h.tot = h1.tot - ( h2.motor ) + h2.tot + fairing.L; %[m] total height
        fairing.h0 = h1.tot - ( h2.motor ) + h2.tot; %[m] 
        h.CG = ( h1.CG.tot * (M1.tot - M1.fairing) + h2.CG.tot * (M2.tot - M2.fairing) +...
                + ( fairing.h0 + (1/3) * fairing.L ) * (fairing.m + M.pay_effective) ) / M.M0;
        h2.h0 = h1.attach - h2.motor; %[m] updated starting height of the second stage
        fairing.h0 = h2.attach; % %[m] updated starting height of the fairing
    
    
        %recover errors:
        err = norm( eps_real(:,i) - [M1.eps; M2.eps] ) / norm([M1.eps; M2.eps]);% + 1e-7 * norm( diam_fairing(i) - diam_fairing(i-1) );
        % err2 = abs( h1.attach - h2.h0 + h2.motor );
        % err = norm( [err1, err2 * 1e-6] );
    
        %update i:
        i = i+1;
        end
        
        if skipcycle%isa(M, "double")
        continue
    end
        
        eps_real(:, i:end) = [];
        eps_end = eps_real(:, end);
        
        %recover loads
        loads.F_drag_tot = sum( loads.F_drag ); %total aerodynamic drag [N]
        
        %mass related parameters
        % M.M0 = M.pay + M1.tot + M2.tot;%[kg] initial mass
        M.M0end = M.M0 - M1.prop;
        % M.M1 = M.M0 - M1.tot; %[kg] mass after first stage separation
        M.M1end = M.M1 - M2.prop;
        M.mr1 = M.M0 / M.M0end; %[-] first stage mass ratio
        M.mr2 = M.M1 / M.M1end; %[-] second stage mass ratio
        M.str1 = M1.str; %[kg] first stage structural mass
        M.str2 = M2.str; %[kg] second stage structural mass
        M.fairing = fairing.m; %[kg] fairing mass
        % M.adapter = adapter.m; %[kg] payload adapter mass
        M.avionics = M1.avionics + M2.avionics; %[kg] total avionics mass
        M.TW1 = T1 * M1.n_mot / ( M.M0 * 9.81 ); %[-] T/W of first stack
        M.TW2 = T2 * M2.n_mot / ( M.M1 * 9.81 ); %[-] T/W of second stack
        M.dv = dv;
        M.err = err;
        M.n_mot1 = M1.n_mot;
        M.n_mot2 = M2.n_mot;
        M.diam1 = diam1;
        M.diam2 = diam2;
        M.th1 = th1;
        M.th2 = th2;
        M.i = i;
        M.stg1 = M1;
        M.stg2 = M2;
        M.wing = M1.wing;
        
        %height
        h.fairing = fairing.L; %[m] height of the fairing
        h.finesse_ratio = h.tot / diam1; %[-] finesse ratio
        h.stg1 = h1;
        h.stg2 = h2;
        
        M_it(m*(j-1)+k) = M;
        h_it(m*(j-1)+k) = h;
        loads_it(m*(j-1)+k) = loads;

    end
end

%plot the two stages:
figure(1);
% fairing.h0 = h.tot - h.fairing;%[m] updated starting height of the fairing
%[~] = fairing_fun1(M.pay, M.pay, fairing, 1);
[~] = fairing_fun1(M.pay, fairing, loads_f, 1);
figure(2);
[~, ~, ~] = inert_mass_common_dome(M2, h2, diam2, AR, loads2, mat2, press2, 1);
figure(3);
[~, ~, ~] = inert_mass_common_dome(M1, h1, diam1, AR, loads1, mat1, press1, 1);

x = [0, diam1/2, diam1/2, diam2/2, diam2/2, 0]';
y = [0, 0, h1.til_tank-h1.dome_rp1, h1.attach, h.tot-2*diam2, h.tot]';
xCG = 0;
yCG = h.CG;
figure(4)
plot([x; -flip(x)], [y; flip(y)], '-k'); grid on; axis equal; hold on;
plot([diam1/2, -diam1/2], [h1.til_tank-h1.dome_rp1,h1.til_tank-h1.dome_rp1], '--k');
plot([diam2/2, -diam2/2], [h1.attach,h1.attach], '--k');
plot([diam2/2, -diam2/2], [h.tot-2*diam2,h.tot-2*diam2], '--k');
plot(xCG, yCG, '+r'); 

save('M_it.mat','M_it');
save('h_it.mat','h_it');

%% Functions

function [m_stag, m_tot, m_prop] = tandem_opt_staging(Is, e, dv, m_pay, fzeroOut)
%This function computes the optimal mass distribution and values between
%stages for a tandem configuration.
% INPUTS: 
% Is : [nx1] [s] vector of the impulses of different stages
% e  : [nx1] [1] vector of the structural mass indexes of different stages
% dv : [1x1] [km/s] target delta_v
% m_pay : [1x1] [kg] payload mass
%
% OUTPUT:
% m_stag : [nx1] [kg] vector of the stages total masses
% m_tot : [1x1] [kg] total initial mass
% m_prop : [nx1] [kg] vector of the stages propellant masses

%get data
g = 9.80665; %[m/s^2]
c = Is*g/1000; %[m/s]

n = length(c);

%set root-finding problem
fun = @(k) dv - c' * log( (c-k) ./ (c.*e) );
lim = min( c-1e-2 );   %min( c.*(1-e));
if fzeroOut == 0
    options = optimset('Display','off');
else
    options = optimset('Display','on');
end

if fun(0) * fun(lim) < 0
    K = fzero(fun, [0, lim], options);
else
    m_stag = 1e12 * ones(n, 1);
    m_tot = sum(m_stag) + m_pay;
    m_str = e .* m_stag;
    m_prop = m_stag - m_str;
    warning('Optimal staging was not possible for this step')
    return
end
%get the payload ratios of each stage ( lam_i = mPL_i / mST_i )
lam = K * ( e ./ ( (1-e).*c - K ) );

%initialize
m_stag = zeros(n, 1);

if n == 1
    m_stag = m_pay / lam;
else
    i = n;
    while i >= 1
        m_stag(i) = ( m_pay + sum(m_stag) ) / lam(i) ;
        i = i-1;
    end
end
%get parameters
m_str = e .* m_stag;
m_prop = m_stag - m_str;
m_tot = sum(m_stag) + m_pay;

end

% function [M, h, th] = inert_mass_common_dome(M, h, diam, AR, loads, mat_id, pressure_type, plotcase)
% 
% % considers thickness equal along the whole tank.
% % safety factor to be defined.
% % the volume is the volume of propellant to be contained: you cannot use
% % this function to evaluate blowdown architectures.
% 
% %constants:
% g = 9.81; %[m/s^2]
% 
% %recover loads:
% nx = loads.nx;
% long_acc = nx*g; %[m/s^2] longitudinal acceleration
% nz = loads.nz;
% tran_acc = nz*g; %[m/s^2] transversal acceleration
% m_sust = loads.m; %[kg] sustained mass
% %h_m_sust = loads.h_m_sust; %[m] height of the lumped mass representing the sustained mass
% 
% %propellant masses
% OF = M.OF;%[-] Ox/Fu ratio
% mlox = M.prop * OF / (1+OF);%[kg] mass of lox
% mrp1 = M.prop * 1  / (1+OF);%[kg] mass of rp1
% 
% %burning time
% M.tb = (mlox+mrp1) / (M.m_dot * M.n_mot); %[s] stage burning time
% 
% %propellant densities
% rholox = M.rholox; %[kg/m^3]
% rhorp1 = M.rhorp1; %[kg/m^3]
% 
% %propellant volumes
% vlox = 1.10 * mlox / rholox; %[m^3] %added 10% margin
% vrp1 = 1.05 * mrp1 / rhorp1; %[m^3] %added 5% margin
% 
% %recover dimensions:
% h_motor = h.motor; %[m] height of the motor
% h0 = h.h0; %[m] height of the bottom part of the stage with the launcher vertically placed on a launch pad 
% R_next = M.R_next; %[m] radius of the subsequent stage
% 
% %get material properties
% mat = mat_switch(mat_id);
% rho = mat.rho;%[kg/m^3] material density
% t_min = mat.t_min; %[m] material manufacturability minimum thickness
% E = mat.E; %[Pa] Young modulus
% sy = mat.sy; %[Pa] yelding stress
% su = mat.su; %[Pa] ultimate stress
% nu = mat.nu; %[-] Poisson's ratio
% 
% %get MEOP (maximum expected operating pressure)
% switch pressure_type 
%     case 0 % case for unpressurized vessel
%         MEOP = 0;
%     case 1 %case for pressure-fed
%         MEOP = 2.8 * 1e6; %[Pa] internal tank pressure
%     case 2 %case for pump-fed
%         MEOP = 500 * 1e3; %[Pa] internal tank pressure
%     case 3 % case for blowdown
%         MEOP = 50 * 1e6; %[Pa] internal tank pressure
% end
% 
% %correction factor
% Km = 1.1; 
% Kp = 1.0;
% MDP = MEOP * Km * Kp;
% jproof = 1.25;
% jburst = 1.5;
% p = MDP * jburst;
% 
% %tanks shape definition
% e = sqrt( 1 - 1 / AR^2 );         %eccentricity of oblate part [-]
% R_int = diam / 2;    %internal radius of the tank [m]
% V_obl = (4/3)*pi*R_int^3 / AR; %volume of the two oblate parts [m^3]
% 
% %check on tank shape
% v_min = min(vlox, vrp1);
% 
% if v_min < (4/3) * pi * R_int^3 / AR  %adopt spherical tanks
%     %lox
%     R_lox = ( 3 * vlox / ( 4*pi ) )^(1/3); %[m] radius of lox spherical tank
%     h_dome_lox = R_lox; %dome height [m]
%     h_cyl_lox = 0;%height of cylindrical part [m]
%     S_lox = 4 * pi * R_lox^2;%surface of tank [m^2]
%     y_lox = 2*R_lox;%fluid level inside the tank [m]
%     AR_lox = 1;
%     h0_lox = (2/3)*R_lox + h0 + h_motor;%[m]
% 
%     %rp1
%     R_rp1 = ( 3 * vrp1 / ( 4*pi ) )^(1/3); %[m] radius of rp1 spherical tank
%     h_dome_rp1 = R_rp1; %dome height [m]
%     h_cyl_rp1 = 0;%height of cylindrical part [m]
%     S_rp1 = 4 * pi * R_rp1^2;%surface of tank [m^2]
%     y_rp1 = 2*R_rp1;%fluid level inside the tank [m]
%     AR_rp1 = 1;
%     h0_rp1 = h0_lox + 2*h_dome_lox + (1/2)*R_rp1;%[m]
%     h0_C1 = h0_rp1 + R_rp1;%[m]
%     h0_C2 = h0_lox + R_lox;%[m]
%     h0_C3 = h0 + h_motor;%[m]
% else  %adopt cyl tanks with shared dome
%     %lox
%     R_lox = R_int; %[m] radius of the lox tank
%     h_dome_lox = R_lox/AR; %dome height [m]
%     V_cyl_lox = vlox - V_obl; %volume of the cylindrical part [m^3]
%     h_cyl_lox = V_cyl_lox / (pi*R_lox^2); %height of cylindrical part [m]
%     y_lox = h_cyl_lox + 2*R_lox/AR;%fluid level inside the tank [m]
%     S_cyl_lox = 2*pi*R_lox* h_cyl_lox; %surface of cylindrical part [m^2]
%     S_obl_lox = 2*pi * R_lox^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
%     S_lox = S_obl_lox + S_cyl_lox; %surface of the lox tank [m^2]
%     AR_lox = AR;
%     h0_lox = (2/3)*R_lox + h0 + h_motor;%[m] 
% 
%     %rp1
%     R_rp1 = R_int; %[m] radius of the rp1 tank
%     h_dome_rp1 = R_rp1/AR; %dome height [m]
%     V_cyl_rp1 = vrp1; %volume of the cylindrical part [m^3]
%     h_cyl_rp1 = V_cyl_rp1 / (pi*R_rp1^2); %height of cylindrical part [m]
%     y_rp1 = h_cyl_rp1 + 2*R_rp1/AR;%fluid level inside the tank [m]
%     S_cyl_rp1 = 2*pi*R_rp1* h_cyl_rp1; %surface of cylindrical part [m^2]
%     S_obl_rp1 = 2*pi * R_rp1^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
%     S_rp1 = S_obl_rp1 + S_cyl_rp1; %surface of the rp1 tank [m^2]
%     AR_rp1 = AR;
%     h0_rp1 = h0_lox + h_dome_lox + h_cyl_lox + 0.16; %[m]. 0.16m added based on Edberg-Costa to account for insulation
%     h0_C1 = h0_rp1 + h_cyl_rp1;%[m]
%     h0_C3 = h0 + h_motor;%[m]
% end
% 
% %recover dimensions
% l_lox = @(t) y_lox + 2*t;%height of the tank [m]
% l_rp1 = @(t) y_rp1 + 2*t;%height of the tank [m]
% h.CG.lox = h_cyl_lox / 2 + h_dome_lox + h0_lox; %[m] COG of lox tank
% h.CG.rp1 = h_cyl_rp1 / 2 + h_dome_rp1 + h0_rp1; %[m] COG of rp1 tank
% M.R_end = R_lox; %[m] ending radius of the stage
% M.R_sta = R_rp1; %[m] starting radius of the stage (w/o considering the top connector) 
% M.diam_cyl = 2 * ( 3*(v_min*AR) / (4*pi) )^(1/3) - 0.01; 
% 
% %pressure at base with longitudinal acceleration (Stevin's law)
% p_lox = p + y_lox * rholox * long_acc; %[Pa] pressure at bottom of tank during acceleration
% p_rp1 = p + y_rp1 * rhorp1 * long_acc; %[Pa] pressure at bottom of tank during acceleration
% 
% %thickness of tanks for pressure
% t_p_lox = diam*p_lox/( 2*sy ); %[m] lox tank thickness
% t_p_rp1 = diam*p_rp1/( 2*sy ); %[m] rp1 tank thickness
% 
% 
% %MASSES ESTIMATION
% 
% %top interstage (first connector) (between top part of the stage and subsequent stage)
% % if R_next > R_rp1  %that is, if  next stage radius is NOT less than this one's
% %     shape1.r = R_rp1; %[m] radius of the cylindrical connector
% %     shape1.h = shape1.r * 5 / 2;
% % else %that is, if  next stage radius is less than this one's
%     shape1.r = [R_next, R_rp1]; %for simplicity we take the same dimensions of the "both-spherical" case
%     shape1.h = shape1.r(2) * 5 / 2;
% % end 
% load1 = loads;
% load1.p = 0; %internal pressure [Pa]
% if nargin < 8
%     [C1.m, C1.th] = buckling_bending(shape1, load1, mat);
% else
%     shape1.h0 = h0_C1;
%     shape1.AR = AR_rp1;
%     [C1.m, C1.th, C1.XY] = buckling_bending(shape1, load1, mat, 1);
%     plot(C1.XY(1,:), C1.XY(2,:), '--k', DisplayName='true');grid on; axis equal; hold on;
% end
% h.CG.C1 = h0_C1 + shape1.h / 2;%[m] exact solution for cyl, approx for cones
% th.C1 = C1.th; %[m] thickness of first connector (interstage)
% 
% 
% %mass estimation of the first tank
% if h_cyl_rp1 > 0    %cyl tank
%     %validate rp1 tank size
%     shape_rp1.r = R_rp1;
%     shape_rp1.h = h_cyl_rp1;
%     load_rp1 = loads;
%     load_rp1.m = m_sust + C1.m + mrp1 * 0.11; %accounts for sustained masses : upper stages, first connector, fairing and rp1 tank structural mass (approximated as 11% of rp1 mass)
%     % load_rp1.F_drag = loads.F_drag; %aerodynamic drag force [N]
%     load_rp1.p = MEOP;%[Pa]
%     [~, t_bb_rp1] = buckling_bending(shape_rp1, load_rp1, mat);
%     t_rp1 = max( t_p_rp1, t_bb_rp1 ); %correction of rp1 tank thickness in case the previous can't sustain the load
% else %spherical tank
%     t_rp1 = max( t_p_rp1, t_min );
% end
% th.rp1 = t_rp1;%[m]
% M.tank_rp1 = rho * S_rp1 * t_rp1; %mass of the empty rp1 tank [kg]
% M.tot_rp1 = M.tank_rp1 + mrp1; %[kg] mass of full rp1 tank
% 
% %middle connector (intertank) (between two tanks)
% if R_rp1 == R_lox %that is, if tanks have the same diameter
%     C2.m = 0;%[kg] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     shape2.h = 0.16; %[m]
%     h.CG.C2 = 0; %[m]
%     th.C2 = 0; %[m] thickness of second connector (intertank)
% else %that is, if at least one tank is spherical
%     shape2.r = [R_rp1, R_lox]; %for simplicity we take the same dimensions of the "both-spherical" case
%     shape2.h = 0.5 * R_rp1 + h_dome_lox + h_dome_rp1;%[m] lenght of the connector between the two tanks
%     load2 = loads;
%     load2.m = m_sust + C1.m + M.tot_rp1; %accounts for sustained masses : upper stages (m_sust), fairing, first connector and rp1 tank mass
%     % load2.F_drag = loads.F_drag; %aerodynamic drag force [N]
%     load2.p = 0; %internal pressure [Pa]
%     if nargin < 8
%         [C2.m, C2.th] = buckling_bending(shape2, load2, mat);
%     else
%         shape2.h0 = h0_C2;
%         shape2.AR = AR_lox;
%         [C2.m, C2.th, C2.XY] = buckling_bending(shape2, load2, mat, 1);
%         plot(C2.XY(1,:), C2.XY(2,:), '--k', DisplayName='true');
%     end
%     h.CG.C2 = h0_C2 + shape2.h / 2;%[m] exact solution for cyl, approx for cones
%     th.C2 = C2.th; %[m] thickness of second connector (intertank)
% end
% 
% %mass estimation of the second tank
% if h_cyl_lox > 0 %cyl tank
%     %validate lox tank size
%     shape_lox.r = R_int;
%     shape_lox.h = h_cyl_lox;
%     load_lox = loads;
%     load_lox.m = m_sust + C1.m + M.tot_rp1 + C2.m + mlox * 0.11; %accounts for sustained masses : upper stages, first and second connector, fairing, rp1 tank and lox tank structural mass (approximated as 11% of lox mass)
%     % load_lox.F_drag = loads.F_drag; %aerodynamic drag force [N]
%     load_lox.p = MEOP; %[Pa] internal pressure
%     [~,t_bb_lox] = buckling_bending(shape_lox, load_lox, mat);
%     t_lox = max( t_p_lox, t_bb_lox ); %correction of lox tank thickness in case the previous can't sustain the load
% else %spherical tank
%     t_lox = max( t_p_lox, t_min);
% end
% th.lox = t_lox;%[m]
% M.lox_insulation = 1.12 * S_lox; %[kg] mass of the insulation layer (from "Launch and Entry Vehicle Design, Univ. Maryland, D.L. Akin")
% M.tank_lox = rho * S_lox * t_lox + M.lox_insulation; %mass of the empty lox tank [kg]
% M.tot_lox = M.tank_lox + mlox; %[kg] mass of full lox tank
% 
% %last connector (aft skirt) (between second tank and motors)
% shape3.r = R_lox;
% shape3.h = (2/3) * R_lox + h_dome_lox; 
% load3 = loads;
% load3.m = m_sust + C1.m + M.tot_rp1 + C2.m + M.tot_lox; %accounts for sustained masses : upper stages, first and second connector, fairing, rp1 tank and lox tank structural mass (approximated as 11% of lox mass)
% % load3.F_drag = loads.F_drag; %aerodynamic force [N]
% load3.p = 0; %internal pressure [Pa]
% if nargin < 8
%     [C3.m, C3.th] = buckling_bending(shape3, load3, mat);
% else
%     shape3.h0 = h0_C3;
%     shape3.AR = AR_lox;
%     [C3.m, C3.th, C3.XY] = buckling_bending(shape3, load3, mat, 1);
%     plot(C3.XY(1,:), C3.XY(2,:), '--k', DisplayName='true');
% end
% C3.m = 0;
% h.CG.C3 = h0_C3 + shape3.h / 2;%[m]
% th.C3 = C3.th; %[m] thickness of third connector (aft skirt)
% 
% %TOTAL MASSES:
% M.tanks = M.tank_lox + M.tank_rp1; %[kg] mass of the two empty tanks
% M.str = M.tanks + M.motor + M.fairing + M.avionics + M.T_struct + M.wiring + M.other + C1.m + C2.m + C3.m; %[kg] inert mass of stage (motors, tanks, connectors, stage fairing, wiring, avionics, structures)
% if M.stg == 1
%     M.recovery = ( 0.07 / ( 1 - 0.07 ) ) * M.str; %[kg]
%     M.str = M.str + M.recovery; %[kg]
% else
%     M.recovery = 0; %[kg]
% end
% M.tot = M.tot_lox + M.tot_rp1 + M.motor + M.fairing + M.avionics + M.T_struct + M.wiring + M.recovery + M.other + C1.m + C2.m + C3.m; %[kg] total mass of motors, tanks, propellant, connectors, wiring, avionics, structures
% M.eps = M.str / M.tot;
% M.C1 = C1;
% M.C2 = C2;
% M.C3 = C3;
% M.mr = ( M.tot + m_sust ) / ( M.str + m_sust ); %[-] approx mass ratio of the stage
% M.n_final = M.Thrust / ( M.str + m_sust ) * 1/g; %[-] ending longitudinal acceleration
% 
% 
% %HEIGHTS
% h.tank_lox = l_lox(t_lox); %[m] height of tank
% h.tank_rp1 = l_rp1(t_rp1); %[m] height of tank
% h.dome_lox = h_dome_lox; %[m]
% h.dome_rp1 = h_dome_rp1; %[m]
% h.R_lox = R_lox; %[m]
% h.R_rp1 = R_rp1; %[m]
% h.C1 = shape1.h; %[m] first connector / top interstage
% h.C2 = shape2.h; %[m] second connector / intertank
% h.C3 = shape3.h; %[m] thirk connector / aft skirt
% h.til_tank = h_motor + h.C3 + h_cyl_lox + h.C2 + h_cyl_rp1 + h_dome_rp1; %[m] height of stage until last tank
% h.tot =      h_motor + h.C3 + h_cyl_lox + h.C2 + h_cyl_rp1 + (5/2)*R_rp1; %[m] total height of stage
% h.CG.avionics = h0_C1 + R_rp1; %[m]
% h.CG.T_struct = h_motor + 0.5 * R_lox; %[m]
% h.CG.tot = (h.CG.lox * M.tot_lox + h.CG.rp1 * M.tot_rp1 +...
%         + h.CG.C1 * C1.m + h.CG.C2 * C2.m + h.CG.C3 * C3.m +...
%         + h.CG.avionics * M.avionics + 0.5 * h_motor * M.motor +...
%         + h.CG.T_struct * M.T_struct) / M.tot;
% h.attach = h0 + h.tot; %[m] height at which the subsequent stage is attached
% 
% %PLOT OF STAGE:
% if nargin > 7
%     if R_lox == R_rp1
%         c = 1;
%     else
%         c = -1;
%     end
%         %lox
%         bottom = @(k) h0_lox + h_dome_lox -sqrt(R_lox^2 - k.^2)/AR_lox;
%         top = @(k) bottom(R_lox) + h_cyl_lox + sqrt(R_lox^2 - k.^2)/AR_lox;
%         K = linspace(-R_lox, R_lox, 1e4);
%         plot(K, bottom(K), 'k'); grid on; axis equal; hold on;
%         plot(K, top(K), 'k');
%         plot([ R_lox,  R_lox], [bottom(R_lox), top(R_lox)], 'k');
%         plot([-R_lox, -R_lox], [bottom(R_lox), top(R_lox)], 'k');
% 
%         %rp1 
%         bottom = @(k) h0_rp1 + c* sqrt(R_rp1^2 - k.^2)/AR_rp1 + 0.5*(1-c)*h_dome_rp1;
%         top = @(k) bottom(R_rp1) + sqrt(R_rp1^2 - k.^2)/AR_rp1 + h_cyl_rp1;
%         K = linspace(-R_rp1, R_rp1, 1e4);
%         plot(K, bottom(K), 'k');
%         plot(K, top(K), 'k');
%         plot([ R_rp1,  R_rp1], [bottom(R_rp1), top(R_rp1)], 'k');
%         plot([-R_rp1, -R_rp1], [bottom(R_rp1), top(R_rp1)], 'k');
% 
%     if h_cyl_lox == 0 %sphere
%         plot(0, h0_lox + h_dome_lox, '+k');
%     end
%     if h_cyl_rp1 == 0 %sphere
%         plot(0, h0_rp1 + h_dome_rp1, '+k');
%     end
%     xlabel('x [m]', 'Interpreter','latex');
%     ylabel('y [m]', 'Interpreter','latex');
% end
% end

% function [m, th, XY] = buckling_bending(shape, load, mat, plotcase)
% 
% % based on NASA papers in shared folder (SP-8007 & SP-8019)
% % computes connectors masses, heights and thicknesses to sustain
% % compression loads and avoid buckling effect
% 
% %constants:
% g = 9.81; %[m/s^2] gravitational acceleration
% 
% %recover loads:
% m = load.m; %sustained mass [kg]
% nx = load.nx; %longitudinal load factor [-]
% K = load.K; %factor of safety [-]
% F_aero = load.F_drag; %aerodynamic drag force [N]
% p = load.p; %internal pressure [Pa]
% M_exp = load.M_exp; %expected flessional 
% 
% %recover material characteristics:
% id = mat.ID; %[-] ID of the material: 1 for Ti, 2 for Al 2XXX, 3 for Steel, 4 for Carbon Fiber, 5 for Al 7XXX
% E = mat.E; %[Pa] Young modulus
% rho = mat.rho; %[kg/m^3]
% t_min1 = mat.t_min; %[m] minimum thickness for manufacturability 
% sy = mat.sy; %[Pa] tensile yield stress
% su = mat.su; %[Pa] tensile ultimate stress
% nu = mat.nu; %[-] Poisson's ratio
% 
% %compute sustained load:
% F_load = m * nx * g + F_aero; %load [N]
% 
% %recover dimensions:
% r = shape.r;
% h = shape.h; %distance between the base of the domes
% if length(r) > 1 %trucated-cone
%     %recover cone shape characteristics
%     alpha = asin( abs( r(2) - r(1) ) / h );
%     L = sqrt( h^2 - (r(2) - r(1))^2 );
%     l = cos(alpha) * L; %height of the shell
%     rho1 = r(1);
%     rho2 = r(2);
%     r(2) = cos(alpha) * r(2); 
%     r(1) = cos(alpha) * r(1);
%     t_min2 = r(2) / 1500; %minimum value for the NASA papers study
%     t_min3a = K * ( r(1) * abs( F_load - p * pi * r(1)^2 ) + 2 * abs( M_exp ) ) /...
%     ( 2 * pi * r(1)^2 * sy ); %minimum value to stay in elastic field of the material (bigger base)
%     t_min3b = K * ( r(2) * abs( F_load - p * pi * r(2)^2 ) + 2 * abs( M_exp ) ) /...
%     ( 2 * pi * r(2)^2 * sy ); %minimum value to stay in elastic field of the material (larger base)
%     t_min3 = min([t_min3a, t_min3b]);%minimum value to stay in elastic field of the material
%     %plot
%     if nargin > 3
%         h0 = shape.h0; %height at which the connector is placed 
%         AR = shape.AR;
%         %for the plotting
%         y = h0 + [rho2*sin(alpha/(2*AR)), rho2*sin(alpha/(2*AR)), rho2*sin(alpha/(2*AR))+l, rho2*sin(alpha/(2*AR))+l, rho2*sin(alpha/(2*AR)), rho2*sin(alpha/(2*AR))];
%         XY = [0, rho2, rho1, -rho1, -rho2, 0; y];
%     end
% else
%     %recover cylinder shape characteristics
%     alpha = 0;
%     L = h;
%     t_min2 = r / 1500; %minimum value for the NASA papers study
%     if nargin > 3
%         h0 = shape.h0; %height at which the connector is placed 
%         %for the plotting
%         y = h0 + [0, 0, L, L, 0, 0];
%         XY = [0, r, r, -r, -r, 0; y];
%     end
%     t_min3 = K * ( r(1) * abs( F_load - p * pi * r(1)^2 ) + 2 * abs( M_exp ) ) /...
%     ( 2 * pi * r(1)^2 * sy ); %minimum value to stay in elastic field of the material
% end
% 
% %CHECK FOR BUCKLING IN AXIAL COMPRESSION + BENDING + INTERNAL PRESSURE:
% %get the wall flexural stiffness per unit width:
% D = @(th) E * th.^3 / ( 12 * (1-nu^2) );
% 
% %compute delta_gamma for pressure-increased performances:
% dg = @(th) d_gamma(p, E, r, th, alpha); 
% 
% 
% %6 SITUATIONS: p OR NON p , CONICAL/CYLINDRICAL (BUT THERE AREN'T CONICAL TANKS), ISOTROPIC/ORTHOTROPIC
% switch id
%     case 600 %(for CF, orthotropic expressions)
% 
%     otherwise
%         if alpha == 0 %(cylindrical shape, NASA SP-8007)
%             k1 = 0.8;
% 
%             %compute knockdown factor:
%             phi = @(th) (1/16) * sqrt(r(1)./th);
%             gP = @(th) 1 - 0.901 * ( 1 - exp( -phi(th) ) );
%             gM = @(th) 1 - 0.731 * ( 1 - exp( -phi(th) ) );
% 
%             %in cyl, distinguish between presurrized and unpressurized cases:
%             if p == 0
%                 k2 = 1;
%             else %( p > 0 )
%                 k2 = 0;
%             end
%         else %(alpha > 0) (conical shape, NASA SP-8019)
%             k1 = 0.5;
%             k2 = 0;
% 
%             %compute knockdown factor:
%             gP = @(th) 0.33; 
%             gM = @(th) 0.41;
%         end
% end
% 
% %compute critical loads expressions
% kx = @(th) k_x(nu, L, min(r), th, gP(th));
% Pcr = @(th) k2 * kx(th) * 2 * pi^3 * D(th) * min(r)   / L^2 +...   %if cyl and p=0 (only metals)
%          + (1-k2) * ( 2*pi    * E * th.^2 .* ( gP(th) ./ sqrt( 3 * (1-nu^2) ) + dg(th) ) + p * pi * min(r)^2 ); %in any other case (only metals)
% Mcr = @(th) k2 * kx(th)   *   pi^3 * D(th) * min(r)^2 / L^2 +...   %if cyl and p=0 (only metals)
%          + (1-k2) * pi*min(r) * E * th.^2 .* ( gM(th) ./ sqrt( 3 * (1-nu^2) ) + dg(th) ) + p * pi * min(r)^2 * k1;%in any other case (only metals)
% 
% %relation to be satisfied in combined stress condition (for AXIAL COMPRESSION + BENDING + INTERNAL PRESSURE)
% f = @(th) K*F_load./Pcr(th) + K*M_exp./Mcr(th) - 1; %this must be <0 to save the structure from failure
% 
% %check if one of the already existing lower bounds satisfies the relation
% %to avoid buckling in axial compression + bending + internal pressure
% t_min = max([t_min1, t_min2, t_min3]);%in order: manufacturability, validity of buckling theory, resistance in plastic field
% f_eval = f(t_min);
% if f_eval < 0
%     th = t_min;%[m]
% else
%     if f(t_min) * f(1) < 0 && imag((t_min)) == 0 && imag(f(1)) == 0
%         th = fzero( f, [t_min, 1]);%[m]
%     else
%         th = 10;
%         warning('For this step it is not possible to avoid buckling');
%     end
% end
% 
% %compute the mass
% if length(r) > 1
%     S = pi * L * ( r(2) + r(1) ); %surface of the truncated cone
%     m = S * th * rho; %[kg]
% else
%     m = 2*pi * r * h * th * rho; %[kg]
% end
% end
% 

function [M, h, th] = inert_mass_common_dome(M, h, diam, AR, loads, mat_id, pressure_type, plotcase)

% considers thickness equal along the whole tank.
% safety factor to be defined.
% the volume is th volume of propellant to be contained: you cannot use
% this function to evaluate blowdown architectures.

%constants:
g = 9.81; %[m/s^2]

%recover loads:
nx = loads.nx;
long_acc = nx*g; %[m/s^2] longitudinal acceleration
nz = loads.nz;
tran_acc = nz*g; %[m/s^2] transversal acceleration
m_sust = loads.m; %[kg] sustained mass
loads.p_hydro = 0; %[Pa] set to zero for connectors 
%h_m_sust = loads.h_m_sust; %[m] height of the lumped mass representing the sustained mass

%propellant masses
OF = M.OF;%[-] Ox/Fu ratio
mlox = M.prop * OF / (1+OF);%[kg] mass of lox
mrp1 = M.prop * 1  / (1+OF);%[kg] mass of rp1

%burning time
M.tb = (mlox+mrp1) / (M.m_dot * M.n_mot); %[s] stage burning time

%propellant densities
rholox = M.rholox; %[kg/m^3]
rhorp1 = M.rhorp1; %[kg/m^3]

%propellant volumes
vlox = 1.04 * mlox / rholox; %[m^3] %added 10% margin
vrp1 = 1.01 * mrp1 / rhorp1; %[m^3] %added 5% margin

%recover dimensions:
h_motor = h.motor; %[m] height of the motor
h0 = h.h0; %[m] height of the bottom part of the stage with the launcher vertically placed on a launch pad 
R_next = M.R_next; %[m] radius of the subsequent stage

%get material properties
mat = mat_switch(mat_id);
rho = mat.rho;%[kg/m^3] material density
t_min = mat.t_min; %[m] material manufacturability minimum thickness
E = mat.E; %[Pa] Young modulus
sy = mat.sy; %[Pa] yelding stress
su = mat.su; %[Pa] ultimate stress
nu = mat.nu; %[-] Poisson's ratio

%get MEOP (maximum expected operating pressure)
switch pressure_type 
    case 0 % case for unpressurized vessel
        MEOP = 0;
    case 1 %case for pressure-fed
        MEOP = 2.8 * 1e6; %[Pa] internal tank pressure
    case 2 %case for pump-fed
        MEOP = 500 * 1e3; %[Pa] internal tank pressure
    case 3 % case for blowdown
        MEOP = 50 * 1e6; %[Pa] internal tank pressure
end

%correction factor
Km = 1.1; 
Kp = 1.0;
MDP = MEOP * Km * Kp;
jproof = 1.25;
jburst = 1.5;
p = MDP * jburst;

%tanks shape definition
e = sqrt( 1 - 1 / AR^2 );         %eccentricity of oblate part [-]
R_int = diam / 2;    %internal radius of the tank [m]
V_obl = (4/3)*pi*R_int^3 / AR; %volume of the two oblate parts [m^3]

%check on tank shape
v_min = min(vlox, vrp1);

if v_min < (4/3) * pi * R_int^3 / AR  %adopt spherical tanks
    %lox
    R_lox = ( 3 * vlox / ( 4*pi ) )^(1/3); %[m] radius of lox spherical tank
    h_dome_lox = R_lox; %dome height [m]
    h_cyl_lox = 0;%height of cylindrical part [m]
    S_lox = 4 * pi * R_lox^2;%surface of tank [m^2]
    y_lox = 2*R_lox;%fluid level inside the tank [m]
    AR_lox = 1;
    h0_lox = (2/3)*R_lox + h0 + h_motor;%[m]

    %rp1
    R_rp1 = ( 3 * vrp1 / ( 4*pi ) )^(1/3); %[m] radius of rp1 spherical tank
    h_dome_rp1 = R_rp1; %dome height [m]
    h_cyl_rp1 = 0;%height of cylindrical part [m]
    S_rp1 = 4 * pi * R_rp1^2;%surface of tank [m^2]
    y_rp1 = 2*R_rp1;%fluid level inside the tank [m]
    AR_rp1 = 1;
    h0_rp1 = h0_lox + 2*h_dome_lox + (1/2)*R_rp1;%[m]
    h0_C1 = h0_rp1 + R_rp1;%[m]
    h0_C2 = h0_lox + R_lox;%[m]
    h0_C3 = h0 + h_motor;%[m]
else  %adopt cyl tanks with shared dome
    %lox
    R_lox = R_int; %[m] radius of the lox tank
    h_dome_lox = R_lox/AR; %dome height [m]
    V_cyl_lox = vlox - V_obl; %volume of the cylindrical part [m^3]
    h_cyl_lox = V_cyl_lox / (pi*R_lox^2); %height of cylindrical part [m]
    y_lox = h_cyl_lox + 2*R_lox/AR;%fluid level inside the tank [m]
    S_cyl_lox = 2*pi*R_lox* h_cyl_lox; %surface of cylindrical part [m^2]
    S_obl_lox = 2*pi * R_lox^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
    S_lox = S_obl_lox + S_cyl_lox; %surface of the lox tank [m^2]
    AR_lox = AR;
    h0_lox = (2/3)*R_lox + h0 + h_motor;%[m] 

    %rp1
    R_rp1 = R_int; %[m] radius of the rp1 tank
    h_dome_rp1 = R_rp1/AR; %dome height [m]
    V_cyl_rp1 = vrp1; %volume of the cylindrical part [m^3]
    h_cyl_rp1 = V_cyl_rp1 / (pi*R_rp1^2); %height of cylindrical part [m]
    y_rp1 = h_cyl_rp1 + 2*R_rp1/AR;%fluid level inside the tank [m]
    S_cyl_rp1 = 2*pi*R_rp1* h_cyl_rp1; %surface of cylindrical part [m^2]
    S_obl_rp1 = 2*pi * R_rp1^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
    S_rp1 = S_obl_rp1 + S_cyl_rp1; %surface of the rp1 tank [m^2]
    AR_rp1 = AR;
    h0_rp1 = h0_lox + h_dome_lox + h_cyl_lox + 0.16; %[m]. 0.16m added based on Edberg-Costa to account for insulation
    h0_C1 = h0_rp1 + h_cyl_rp1;%[m]
    h0_C3 = h0 + h_motor;%[m]
end

%recover dimensions
l_lox = @(t) y_lox + 2*t;%height of the tank [m]
l_rp1 = @(t) y_rp1 + 2*t;%height of the tank [m]
h.CG.lox = h_cyl_lox / 2 + h_dome_lox + h0_lox; %[m] COG of lox tank
h.CG.rp1 = h_cyl_rp1 / 2 + h_dome_rp1 + h0_rp1; %[m] COG of rp1 tank
M.R_end = R_lox; %[m] ending radius of the stage
M.R_sta = R_rp1; %[m] starting radius of the stage (w/o considering the top connector) 
diam_cyl_max = 2 * ( 3*(v_min*AR) / (4*pi) )^(1/3) - 0.01;
% if diam_cyl_max > R_int
%     M.diam_cyl = R_int;
% else
    M.diam_cyl = diam_cyl_max;
% end

%pressure at base with longitudinal acceleration (Stevin's law)
p_lox = p + y_lox * rholox * long_acc; %[Pa] pressure at bottom of tank during acceleration
p_rp1 = p + y_rp1 * rhorp1 * long_acc; %[Pa] pressure at bottom of tank during acceleration

%thickness of tanks for pressure
t_p_lox = diam*p_lox/( 2*sy ); %[m] lox tank thickness
t_p_rp1 = diam*p_rp1/( 2*sy ); %[m] rp1 tank thickness

%MASSES ESTIMATION

%top interstage (first connector) (between top part of the stage and subsequent stage)
% if R_next > R_rp1  %that is, if  next stage radius is NOT less than this one's
%     shape1.r = R_rp1; %[m] radius of the cylindrical connector
%     shape1.h = shape1.r * 5 / 2;
% else %that is, if  next stage radius is less than this one's
    % shape1.r = [R_next, R_rp1]; %for simplicity we take the same dimensions of the "both-spherical" case
    % shape1.h = shape1.r(2) * 5 / 2;
% end 
shape1.r = [R_next, R_rp1]; %for simplicity we take the same dimensions of the "both-spherical" case
if M.stg == 1
    shape1.h = h_dome_rp1 + h.motor_next + 0.1;
else
    shape1.h = (2/3)*R_rp1 + h_dome_rp1;
end
load1 = loads;
load1.p = 0; %internal pressure [Pa]
if nargin < 8
    [C1.m, C1.th] = buckling_bending(shape1, load1, mat);
else
    shape1.h0 = h0_C1;
    shape1.AR = AR_rp1;
    [C1.m, C1.th, C1.XY] = buckling_bending(shape1, load1, mat, 1);
    plot(C1.XY(1,:), C1.XY(2,:), '--k', DisplayName='true');grid on; axis equal; hold on;
end
h.CG.C1 = h0_C1 + shape1.h / 2;%[m] exact solution for cyl, approx for cones
th.C1 = C1.th; %[m] thickness of first connector (interstage)


%mass estimation of the first tank
if h_cyl_rp1 > 0    %cyl tank
    %validate rp1 tank size
    shape_rp1.r = R_rp1;
    shape_rp1.h = h_cyl_rp1;
    load_rp1 = loads;
    load_rp1.p_hydro = y_lox * rholox * long_acc; %[Pa] hydrostatic pressure due to liquid
    load_rp1.m = m_sust + C1.m + mrp1 * 0.11; %accounts for sustained masses : upper stages, first connector, fairing and rp1 tank structural mass (approximated as 11% of rp1 mass)
    % load_rp1.F_drag = loads.F_drag; %aerodynamic drag force [N]
    load_rp1.p = MEOP;%[Pa]
    [~, t_bb_rp1] = buckling_bending(shape_rp1, load_rp1, mat);
    t_rp1 = t_bb_rp1; %correction of rp1 tank thickness in case the previous can't sustain the load
else %spherical tank
    t_rp1 = max( t_p_rp1, t_min );
end
th.rp1 = t_rp1;%[m]
M.tank_rp1 = rho * S_rp1 * t_rp1; %mass of the empty rp1 tank [kg]
M.tot_rp1 = M.tank_rp1 + mrp1; %[kg] mass of full rp1 tank

%middle connector (intertank) (between two tanks)
if R_rp1 == R_lox %that is, if tanks have the same diameter
    C2.m = 0;%[kg] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    shape2.h = 0.16; %[m]
    h.CG.C2 = 0; %[m]
    th.C2 = 0; %[m] thickness of second connector (intertank)
else %that is, if at least one tank is spherical
    shape2.r = [R_rp1, R_lox]; %for simplicity we take the same dimensions of the "both-spherical" case
    shape2.h = 0.5 * R_rp1 + h_dome_lox + h_dome_rp1;%[m] lenght of the connector between the two tanks
    load2 = loads;
    load2.m = m_sust + C1.m + M.tot_rp1; %accounts for sustained masses : upper stages (m_sust), fairing, first connector and rp1 tank mass
    % load2.F_drag = loads.F_drag; %aerodynamic drag force [N]
    load2.p = 0; %internal pressure [Pa]
    if nargin < 8
        [C2.m, C2.th] = buckling_bending(shape2, load2, mat);
    else
        shape2.h0 = h0_C2;
        shape2.AR = AR_lox;
        [C2.m, C2.th, C2.XY] = buckling_bending(shape2, load2, mat, 1);
        plot(C2.XY(1,:), C2.XY(2,:), '--k', DisplayName='true');
    end
    h.CG.C2 = h0_C2 + shape2.h / 2;%[m] exact solution for cyl, approx for cones
    th.C2 = C2.th; %[m] thickness of second connector (intertank)
end

%mass estimation of the second tank
if h_cyl_lox > 0 %cyl tank
    %validate lox tank size
    shape_lox.r = R_int;
    shape_lox.h = h_cyl_lox;
    load_lox = loads;
    load_lox.p_hydro = y_rp1 * rhorp1 * long_acc; %[Pa] hydrostatic pressure due to liquid
    load_lox.m = m_sust + C1.m + M.tot_rp1 + C2.m + mlox * 0.11; %accounts for sustained masses : upper stages, first and second connector, fairing, rp1 tank and lox tank structural mass (approximated as 11% of lox mass)
    % load_lox.F_drag = loads.F_drag; %aerodynamic drag force [N]
    load_lox.p = MEOP; %[Pa] internal pressure
    [~,t_bb_lox] = buckling_bending(shape_lox, load_lox, mat);
    t_lox = t_bb_lox; %correction of lox tank thickness in case the previous can't sustain the load
else %spherical tank
    t_lox = max( t_p_lox, t_min);
end
th.lox = t_lox;%[m]
M.lox_insulation = 1.12 * S_lox; %[kg] mass of the insulation layer (from "Launch and Entry Vehicle Design, Univ. Maryland, D.L. Akin")
M.tank_lox = rho * S_lox * t_lox + M.lox_insulation; %mass of the empty lox tank [kg]
M.tot_lox = M.tank_lox + mlox; %[kg] mass of full lox tank

%last connector (aft skirt) (between second tank and motors)
% shape3.r = R_lox;
% shape3.h = (2/3) * R_lox + h_dome_lox; 
shape3.r = R_lox;
if M.stg == 1
    shape3.h = 2*R_lox;
else
    shape3.h = (2/3)*R_lox + h_dome_lox;
end
load3 = loads;
load3.m = m_sust + C1.m + M.tot_rp1 + C2.m + M.tot_lox; %accounts for sustained masses : upper stages, first and second connector, fairing, rp1 tank and lox tank structural mass (approximated as 11% of lox mass)
% load3.F_drag = loads.F_drag; %aerodynamic force [N]
load3.p = 0; %internal pressure [Pa]
if nargin < 8
    [C3.m, C3.th] = buckling_bending(shape3, load3, mat);
else
    shape3.h0 = h0_C3;
    shape3.AR = AR_lox;
    [C3.m, C3.th, C3.XY] = buckling_bending(shape3, load3, mat, 1);
    plot(C3.XY(1,:), C3.XY(2,:), '--k', DisplayName='true');
end
%C3.m = 0;
h.CG.C3 = h0_C3 + shape3.h / 2;%[m]
th.C3 = C3.th; %[m] thickness of third connector (aft skirt)

%TOTAL MASSES:
M.tanks = M.tank_lox + M.tank_rp1; %[kg] mass of the two empty tanks
M.str = M.tanks + M.motor + M.fairing + M.avionics + M.T_struct + M.wiring + M.other + M.wing + C1.m + C2.m + C3.m; %[kg] inert mass of stage (motors, tanks, connectors, stage fairing, wiring, avionics, structures)
if M.stg == 1
    M.recovery = ( 0.07 / ( 1 - 0.07 ) ) * M.str; %[kg]
    M.str = M.str + M.recovery; %[kg]
else
    M.recovery = 0; %[kg]
end
M.tot = M.tot_lox + M.tot_rp1 + M.motor + M.fairing + M.avionics + M.T_struct + M.wiring + M.recovery + M.other + M.wing + C1.m + C2.m + C3.m; %[kg] total mass of motors, tanks, propellant, connectors, wiring, avionics, structures
M.eps = M.str / M.tot;
M.C1 = C1;
M.C2 = C2;
M.C3 = C3;
M.skirts = C1.m + C2.m + C3.m; 
M.tank_no_ins = M.tanks - M.lox_insulation; 
M.mr = ( M.tot + m_sust ) / ( M.str + m_sust ); %[-] approx mass ratio of the stage
M.n_final = M.Thrust / ( M.str + m_sust ) * 1/g; %[-] ending longitudinal acceleration


%HEIGHTS
h.tank_lox = l_lox(t_lox); %[m] height of tank
h.tank_rp1 = l_rp1(t_rp1); %[m] height of tank
h.dome_lox = h_dome_lox; %[m]
h.dome_rp1 = h_dome_rp1; %[m]
h.R_lox = R_lox; %[m]
h.R_rp1 = R_rp1; %[m]
h.C1 = shape1.h; %[m] first connector / top interstage
h.cyl_rp1 = h_cyl_rp1; %[m] rp1 tank cylindrical part
h.C2 = shape2.h; %[m] second connector / intertank
h.cyl_lox = h_cyl_lox; %[m] lox tank cylindrical part
h.C3 = shape3.h; %[m] thirk connector / aft skirt
h.til_tank = h_motor + h.C3 + h_cyl_lox + h.C2 + h_cyl_rp1 + h_dome_rp1; %[m] height of stage until last tank
h.tot =      h_motor + h.C3 + h_cyl_lox + h.C2 + h_cyl_rp1 + h.C1; %[m] total height of stage
h.CG.avionics = h0_C1 + R_rp1; %[m]
h.CG.T_struct = h_motor + 0.5*(h.C3-h_dome_lox); %[m]
h.CG.tot = (h.CG.lox * M.tot_lox + h.CG.rp1 * M.tot_rp1 +...
        + h.CG.C1 * C1.m + h.CG.C2 * C2.m + h.CG.C3 * C3.m +...
        + h.CG.avionics * M.avionics + 0.5 * h_motor * M.motor +...
        + h.CG.T_struct * M.T_struct) / M.tot;
h.attach = h0 + h.tot; %[m] height at which the subsequent stage is attached

%PLOT OF STAGE:
if nargin > 7
    if R_lox == R_rp1
        c = 1;
    else
        c = -1;
    end
        %lox
        bottom = @(k) h0_lox + h_dome_lox -sqrt(R_lox^2 - k.^2)/AR_lox;
        top = @(k) bottom(R_lox) + h_cyl_lox + sqrt(R_lox^2 - k.^2)/AR_lox;
        K = linspace(-R_lox, R_lox, 1e4);
        plot(K, bottom(K), 'k'); grid on; axis equal; hold on;
        plot(K, top(K), 'k');
        plot([ R_lox,  R_lox], [bottom(R_lox), top(R_lox)], 'k');
        plot([-R_lox, -R_lox], [bottom(R_lox), top(R_lox)], 'k');
        
        %rp1 
        bottom = @(k) h0_rp1 + c* sqrt(R_rp1^2 - k.^2)/AR_rp1 + 0.5*(1-c)*h_dome_rp1;
        top = @(k) bottom(R_rp1) + sqrt(R_rp1^2 - k.^2)/AR_rp1 + h_cyl_rp1;
        K = linspace(-R_rp1, R_rp1, 1e4);
        plot(K, bottom(K), 'k');
        plot(K, top(K), 'k');
        plot([ R_rp1,  R_rp1], [bottom(R_rp1), top(R_rp1)], 'k');
        plot([-R_rp1, -R_rp1], [bottom(R_rp1), top(R_rp1)], 'k');

    if h_cyl_lox == 0 %sphere
        plot(0, h0_lox + h_dome_lox, '+k');
    end
    if h_cyl_rp1 == 0 %sphere
        plot(0, h0_rp1 + h_dome_rp1, '+k');
    end
    xlabel('x [m]', 'Interpreter','latex');
    ylabel('y [m]', 'Interpreter','latex');
end
end

function [m, th, XY] = buckling_bending(shape, load, mat, plotcase)

% based on NASA papers in shared folder (SP-8007 & SP-8019)
% computes connectors masses, heights and thicknesses to sustain
% compression loads and avoid buckling effect

%constants:
g = 9.81; %[m/s^2] gravitational acceleration

%recover loads:
m = load.m; %sustained mass [kg]
nx = load.nx; %longitudinal load factor [-]
K = load.K; %factor of safety [-]
F_aero = load.F_drag; %aerodynamic drag force [N]
p = load.p; %internal pressure [Pa]
p_hydro = load.p_hydro; %hydrostatic pressure [Pa]
M_exp = load.M_exp; %expected flessional 

%recover material characteristics:
id = mat.ID; %[-] ID of the material: 1 for Ti, 2 for Al 2XXX, 3 for Steel, 4 for Carbon Fiber, 5 for Al 7XXX
E = mat.E; %[Pa] Young modulus
rho = mat.rho; %[kg/m^3]
t_min1 = mat.t_min; %[m] minimum thickness for manufacturability 
sy = mat.sy; %[Pa] tensile yield stress
su = mat.su; %[Pa] tensile ultimate stress
nu = mat.nu; %[-] Poisson's ratio

%compute sustained load:
F_load = m * nx * g + F_aero; %load [N]

%recover dimensions:
r = shape.r;
h = shape.h; %distance between the base of the domes
if length(r) > 1 %trucated-cone
    %recover cone shape characteristics
    alpha = asin( abs( r(2) - r(1) ) / h );
    L = sqrt( h^2 - (r(2) - r(1))^2 );
    l = cos(alpha) * L; %height of the shell
    rho1 = r(1);
    rho2 = r(2);
    r(2) = cos(alpha) * r(2); 
    r(1) = cos(alpha) * r(1);
    t_min2 = r(2) / 1500; %minimum value for the NASA papers study
    %check for P, M, p (top of tank)
    t_min3a = K * ( r(1) * abs( F_load - p * pi * r(1)^2 ) + 2 * abs( M_exp ) ) /...
    ( 2 * pi * r(1)^2 * sy ); %minimum value to stay in elastic field of the material (bigger base)
    t_min3b = K * ( r(2) * abs( F_load - p * pi * r(2)^2 ) + 2 * abs( M_exp ) ) /...
    ( 2 * pi * r(2)^2 * sy ); %minimum value to stay in elastic field of the material (larger base)
    t_min3 = min([t_min3a, t_min3b]);%minimum value to stay in elastic field of the material
    %check for P, M, p, p_hydro (bottom of tank)
    t_min4a = K * ( r(1) * abs( F_load - ( p + p_hydro ) * pi * r(1)^2 ) + 2 * abs( M_exp ) ) /...
    ( 2 * pi * r(1)^2 * sy ); %minimum value to stay in elastic field of the material (bigger base)
    t_min4b = K * ( r(2) * abs( F_load - ( p + p_hydro ) * pi * r(2)^2 ) + 2 * abs( M_exp ) ) /...
    ( 2 * pi * r(2)^2 * sy ); %minimum value to stay in elastic field of the material (larger base)
    t_min4 = min([t_min4a, t_min4b]);%minimum value to stay in elastic field of the material
    %check for p, p_hydro (hoop)
    t_min5 = K * abs( p + p_hydro ) * max(r) / sy; %minimum value to stay in elastic field of the material
    %check for M (horizontal unfilled rocket, p=p_hydro=0)
    t_min6 = K * abs( M_exp / 2 ) / ( pi * min(r)^2 * sy ); %minimum value to stay in elastic field of the material
    %plot
    if nargin > 3
        h0 = shape.h0; %height at which the connector is placed 
        AR = shape.AR;
        %for the plotting
        y = h0 + [rho2*sin(alpha/(2*AR)), rho2*sin(alpha/(2*AR)), rho2*sin(alpha/(2*AR))+l, rho2*sin(alpha/(2*AR))+l, rho2*sin(alpha/(2*AR)), rho2*sin(alpha/(2*AR))];
        XY = [0, rho2, rho1, -rho1, -rho2, 0; y];
    end
else
    %recover cylinder shape characteristics
    alpha = 0;
    L = h;
    t_min2 = r / 1500; %minimum value for the NASA papers study
    if nargin > 3
        h0 = shape.h0; %height at which the connector is placed 
        %for the plotting
        y = h0 + [0, 0, L, L, 0, 0];
        XY = [0, r, r, -r, -r, 0; y];
    end
    %check for P, M, p (top of tank)
    t_min3 = K * ( r(1) * abs( F_load - p * pi * r(1)^2 ) + 2 * abs( M_exp ) ) /...
    ( 2 * pi * r(1)^2 * sy ); %minimum value to stay in elastic field of the material
    %check for P, M, p, p_hidro (bottom of tank)
    t_min4 = K * ( r(1) * abs( F_load - ( p + p_hydro ) * pi * r(1)^2 ) + 2 * abs( M_exp ) ) /...
    ( 2 * pi * r(1)^2 * sy ); %minimum value to stay in elastic field of the material
    %check for p, p_hydro (tank hoop resistance)
    t_min5 = K * abs( p + p_hydro ) * r(1) / sy; %minimum value to stay in elastic field of the material
    %check for M (horizontal unfilled rocket, p=p_hydro=0)
    t_min6 = K * abs( M_exp / 2 ) / ( pi * r(1)^2 * sy ); %minimum value to stay in elastic field of the material
end

%CHECK FOR BUCKLING IN AXIAL COMPRESSION + BENDING + INTERNAL PRESSURE:
%get the wall flexural stiffness per unit width:
D = @(th) E * th.^3 / ( 12 * (1-nu^2) );

%compute delta_gamma for pressure-increased performances:
dg = @(th) d_gamma(p, E, r, th, alpha); 


%6 SITUATIONS: p OR NON p , CONICAL/CYLINDRICAL (BUT THERE AREN'T CONICAL TANKS), ISOTROPIC/ORTHOTROPIC
switch id
    case 600 %(for CF, orthotropic expressions)

    otherwise
        if alpha == 0 %(cylindrical shape, NASA SP-8007)
            k1 = 0.8;

            %compute knockdown factor:
            phi = @(th) (1/16) * sqrt(r(1)./th);
            gP = @(th) 1 - 0.901 * ( 1 - exp( -phi(th) ) );
            gM = @(th) 1 - 0.731 * ( 1 - exp( -phi(th) ) );

            %in cyl, distinguish between presurrized and unpressurized cases:
            if p == 0
                k2 = 1;
            else %( p > 0 )
                k2 = 0;
            end
        else %(alpha > 0) (conical shape, NASA SP-8019)
            k1 = 0.5;
            k2 = 0;

            %compute knockdown factor:
            gP = @(th) 0.33; 
            gM = @(th) 0.41;
        end
end

%compute critical loads expressions
kxP = @(th) k_x(nu, L, min(r), th, gP(th));
kxM = @(th) k_x(nu, L, min(r), th, gM(th));
Pcr = @(th) k2 * kxP(th) * 2 * pi^3 * D(th) * min(r)   / L^2 +...   %if cyl and p=0 (only metals)
         + (1-k2) * ( 2*pi    * E * th.^2 .* ( gP(th) ./ sqrt( 3 * (1-nu^2) ) + dg(th) ) + p * pi * min(r)^2 ); %in any other case (only metals)
Mcr = @(th) k2 * kxM(th)   *   pi^3 * D(th) * min(r)^2 / L^2 +...   %if cyl and p=0 (only metals)
         + (1-k2) * pi*min(r) * E * th.^2 .* ( gM(th) ./ sqrt( 3 * (1-nu^2) ) + dg(th) ) + p * pi * min(r)^2 * k1;%in any other case (only metals)

Mcr_unp = @(th) kxM(th) * pi^3 * D(th) * min(r)^2 / L^2;   %if cyl and p=0 (only metals)

%relation to be satisfied in combined stress condition (for AXIAL COMPRESSION + BENDING + INTERNAL PRESSURE)
f1 = @(th) K*F_load./Pcr(th) + K*M_exp./Mcr(th) - 1; %this must be <0 to save the structure from failure

%relation to be satisfied in horizontal, unfilled rocket (ONLY BENDING)
f2 = @(th) K*M_exp - Mcr_unp(th); %this must be <0 to save the structure from failure

%check if one of the already existing lower bounds satisfies the relation
%to avoid buckling in axial compression + bending + internal pressure
t_min = max([t_min1, t_min2, t_min3, t_min4, t_min5, t_min6]);%in order: manufacturability, validity of buckling theory, resistance in plastic field (for P+M+p, P+M+p+p_hydro, hoop, only M) 
f_eval1 = f1(t_min);
f_eval2 = f2(t_min);
if f_eval1 < 0 && f_eval2 < 0
    th = t_min;%[m]
else
    if imag(f_eval1) == 0 && imag(f_eval2) == 0
        if f_eval1 < 0 && f_eval2 > 0
            th = fzero( f2, [t_min, 1e6]);%[m]
        elseif f_eval1 > 0 && f_eval2 < 0
            th = fzero( f1, [t_min, 1e6]);%[m]
        else
            th1 = fzero( f1, [t_min, 1e6]);%[m]
            th2 = fzero( f2, [t_min, 1e6]);%[m]
            th = max( th1, th2 );%[m]
        end
    else
        th = 10;
        warning('For this step it is not possible to avoid buckling');
    end
end

%compute the mass
if length(r) > 1
    S = pi * L * ( r(2) + r(1) ); %surface of the truncated cone
    m = S * th * rho; %[kg]
else
    m = 2*pi * r * h * th * rho; %[kg]
end
end

function dg = d_gamma(p, E, r, th, alpha)

% based on NASA paper SP-8007-2020/REV 2:

if p == 0
    dg = 0;
else
    if nargin < 5
        alpha = 0;
    end
    param = (p/E)*(r(1)/ (th*cos(alpha)) )^2;

    % from curve fitting with 
    % param = [0.02, 0.1, 1, 10]
    % d_gamma = [0.026, 0.09, 0.2, 0.23]
    % and f(x) = (a*x+b)/(c*x+d)
    % we obtained the fitting (R-square = 0.99995, SSE = 1.241e-6, DFE = 0):

    a  = 0.3045;   
    b  = 0.0001;
    c  = 1.3064;
    d  = 0.2104;

    dg = ( a * param + b ) / ( c * param + d );
end
end

function kx = k_x(nu, L, r, th, g)

% based on NASA paper SP-8007-2020/REV 2:

gZ = @(t) g * L^2 / (r*t) * sqrt(1-nu^2);

% from curve fitting with 
% param = [0.1, 1, 3, 100, 1000]
% d_gamma = [1, 1.1, 2, 64, 700]
% and f(x) = c*(a+x^2)^b
% we obtained the fitting (R-square = 1.0000, SSE = 0.0025, DFE = 1):
 
a = 3.2513;    
b = 0.5195;   
c = 0.5348;

kx = c * (a + gZ(th)^2 ) ^ b;
end

% function [fairing] = fairing_fun(m_pay_max, m_pay, fairing, plotcase)
% 
% % computes the shape and mass of the fairing.
% % based on cubesats volumetric density constraints, takes as input also the
% % aerodynamic loads.
% % fairing is simplified as a cylinder and a cone.
% 
% %recover fixed shape characteristics
% d0 = fairing.base_diam; %[m] diameter of the base
% L_nose = fairing.base_diam * 2; %[m] length of the conical nose (from Edberg-Costa)
% 
% %cubesats characteristics:
% vol_den = 1.33 * 1e3; %[kg/m^3] volumetric density of cubesats
% a = 0.1; %[m] cubesat unit edge length 
% b = 0.05; %[m] margin distance between payload and fairing
% ha = 0.40; %[m] adapter assumed height; (suggested by electron pl user guide)
% ra = d0 * 5 / 12; %[m] adapter initial radius
% Ra = ra * 2 / 3; %[m] adapter final radius
% 
% %compute payload volume (considering cubesats sizes)
% V = m_pay_max / vol_den; %[m^3] volume of max payload
% V_real = m_pay / vol_den; %[m^3] volume of real payload
% 
% %find usable base for payload
% y = 0 : a : d0/2; 
% N = length(y) - 1; 
% base = 0;
% x = zeros(N, 1);
% for j = 1 : N
%     x(j) = floor( sqrt( 0.25*d0^2 - y(j)^2 - b ) / a ) * a; %[m] at y(i) level, how much space i have to fill with cubesats
%     base = base + 4 * a * x(j); %[m^2] update base value 
% end
% L_min = V / base; %[m] height of the ammissible payload
% L_real = V_real / base; %[m] height of the real payload 
% 
% %find height of the cylindrical part of the fairing
% L1 = L_min + b + ha; %[m] 
% 
% %surface of the fairing
% S_cyl = L1 * pi * d0; %[m^2] cylinder surface
% S_cone = pi * sqrt( d0^2 / 4 + L_nose^2 ) * d0/2; %[m^2] conical surface
% fairing.S = S_cyl + S_cone; %[m^2] total surface
% 
% %compute mass (from "Launch and Entry Vehicle Design, Univ. Maryland, D.L.Akin")
% if fairing.mat_id == 4 %(composite)
%     fairing.m = 9.89 * fairing.S;%[kg] from edberg-costa
% else
%     fairing.m = 13.3 * fairing.S;%[kg] from edberg-costa
% end
% 
% %retrieve fairing parameters
% % fairing.M = nose.M + cyl.M; %[kg] total mass
% fairing.L = L1 + L_nose; %[m] total height
% fairing.d = d0; %[m] maximum diameter
% % fairing.th_cone = nose.th; %[m] thickness of the nose
% % fairing.th_cyl = cyl.th; %[m] thickness of the cylinder
% fairing.V = V; %[m^3] maximum volume for the payload
% 
% %plot of the fairing
% if nargin > 3
%     h0 = fairing.h0;
%     X = [d0/2, d0/2, 0, -d0/2, -d0/2];
%     Y = h0 + [0, L1, fairing.L, L1, 0];
%     plot(X, Y, 'k', DisplayName='true'); grid on; hold on; axis equal;%fairing
%     X_a = [0, ra, Ra, -Ra, -ra, 0];
%     Y_a = h0 + [0, 0, ha, ha, 0, 0];
%     plot(X_a, Y_a, '--r', DisplayName='true'); %adapter
%     X_p = [0, x(1), x(1), -x(1), -x(1), 0];
%     Y_p = h0 + [ha, ha, ha + L_real, ha + L_real, ha, ha]; 
%     plot(X_p, Y_p, 'b', DisplayName='true'); %payload
%     X_r = [0, 0.6, 0.6, -0.6, -0.6, 0];
%     Y_r = [0, 0, h0, h0, 0, 0];
%     plot(X_r, Y_r, '--k', DisplayName='true');
%     legend('Fairing', 'Adapter', 'Payload','Hypothetical Rocket', 'interpreter', 'latex');
%     xlabel('x [m]', 'Interpreter','latex');
%     ylabel('y [m]', 'Interpreter','latex');
% end
% end

% function [adapter] = adapter_fun(adapter, loads)
% 
% % this function computer adapter characteristics.
% % it uses the loads, the fixed maximum diameter, adn the material
% 
% %recover fixed shape characteristics
% d0 = adapter.base_diam; %[m] diameter of the base
% adapter.r = [d0/3, d0/2]; %[m] payload adapter r1, r2
% adapter.h = d0 / 4; %[m] payload adapter height 
% 
% %update internal pressure load:
% loads.p = 0; %[Pa]
% 
% %recover material properties
% mat = mat_switch(adapter.mat_id);
% 
% %compute adapter characteristics
% [adapter.m , adapter.th] = buckling_bending(adapter, loads, mat);
% 
% end
% 

function [mat] = mat_switch(mat_id)

%this function gets the info about the selected material:
% mat_id = 1 for Ti, 2 for Al 2XXX, 3 for Steel, 4 for Carbon Fiber, 5 for Al 7XXX

%switch
switch mat_id 
    case 1 % Ti6Al4V
        rho = 4500; %[kg/m^3]
        t_min = 0.5 * 1e-3; %[m] minimum thickness for manufacturability
        E = 110 * 1e9; %[Pa] young modulus
        sy = 900 * 1e6; %[Pa] tensile yield stress
        su = 950 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.34; %[-] Poisson's ratio
    case 2 % Al 2XXX
        rho = 2700; %[kg/m^3]
        t_min = 0.5 * 1e-3; %[m] minimum thickness for manufacturability
        E = 70 * 1e9; %[Pa] young modulus
        sy = 290 * 1e6; %[Pa] tensile yield stress
        su = 390 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.33; %[-] Poisson's ratio
    case 3 % Steel
        rho = 7800; %[kg/m^3]
        t_min = 0.25 * 1e-3; %[m] minimum thickness for manufacturability
        E = 200 * 1e9; %[Pa] young modulus
        sy = 350 * 1e6; %[Pa] tensile yield stress
        su = 420 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.27; %[-] Poisson's ratio
    case 4 % Carbon fiber  Toray M46J
        rho = 1600; %[kg/m^3]
        t_min = 0.90 * 1e-3; %[m] minimum thickness for manufacturability
        E = 222 * 1e9; %[Pa] young modulus
        sy = 1090 * 1e6; %[Pa] tensile yield stress
        su = sy; %[Pa] tensile ultimate stress
        nu = 0.28; %[-] Poisson's ratio
    case 5 % Al 7075 T6
        rho = 2810; %[kg/m^3]
        t_min = 1.06 * 1e-3; %[m] minimum thickness for manufacturability
        E = 70 * 1e9; %[Pa] young modulus
        sy = 440 * 1e6; %[Pa] tensile yield stress
        su = 517 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.33; %[-] Poisson's ratio
    case 6 % AlLi (2090)
        rho = 2590; %[kg/m^3]
        t_min = 0.5 * 1e-3; %[m] minimum thickness for manufacturability  
        E = 76 * 1e9; %[Pa] young modulus
        sy = 500 * 1e6; %[Pa] tensile yield stress
        su = 550 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.34; %[-] Poisson's ratio
    case 7 % Carbon fiber Hexcel® HexTow® IM7
        rho = 1600; %[kg/m^3]
        t_min = 0.90 * 1e-3; %[m] minimum thickness for manufacturability
        E = 148 * 1e9; %[Pa] young modulus
        sy = 1000 * 1e6; %[Pa] tensile yield stress
        su = sy; %[Pa] tensile ultimate stress
        nu = 0.28; %[-] Poisson's ratio
    case 8 % Al 6061 T6
        rho = 2700; %[kg/m^3]
        t_min = 0.5 * 1e-3; %[m] minimum thickness for manufacturability
        E = 68.9 * 1e9; %[Pa] young modulus
        sy = 276 * 1e6; %[Pa] tensile yield stress
        su = 310 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.33; %[-] Poisson's ratio
    case 9 % 300M Steel alloy
        rho = 7830; %[kg/m^3]
        t_min = 0.25 * 1e-3; %[m] minimum thickness for manufacturability
        E = 207 * 1e9; %[Pa] young modulus
        sy = 1586 * 1e6; %[Pa] tensile yield stress
        su = 1931 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.3; %[-] Poisson's ratio
    case 10 % 2219 Al alloy
        rho = 2840; %[kg/m^3]
        t_min = 0.5 * 1e-3; %[m] minimum thickness for manufacturability
        E = 73.1 * 1e9; %[Pa] young modulus
        sy = 350 * 1e6; %[Pa] tensile yield stress
        su = 440 * 1e6; %[Pa] tensile ultimate stress
        nu = 0.33; %[-] Poisson's ratio
    case 11 % Carbon fiber  (MatWeb)
        rho = 1420; %[kg/m^3]
        t_min = 3 * 1e-3; %[m] minimum thickness for manufacturability
        E = 101 * 1e9; %[Pa] young modulus
        sy = 1260 * 1e6; %[Pa] tensile yield stress
        su = sy; %[Pa] tensile ultimate stress
        nu = 0.286; %[-] Poisson's ratio
end

%recover material properties:
mat.ID = mat_id; %[-] ID of the material: 1 for Ti, 2 for Al 2XXX, 3 for Steel, 4 for Carbon Fiber, 5 for Al 7XXX
mat.rho = rho; %[kg/m^3] material density
mat.t_min = t_min; %[m] material manufacturability minimum thickness
mat.E = E; %[Pa] Young modulus
mat.sy = sy; %[Pa] yelding stress
mat.su = su; %[Pa] ultimate stress
mat.nu = nu; %[-] Poisson's ratio

end

function [adapter] = adapter_fun1(adapter, loads)

% this function computer adapter characteristics.
% it uses the loads, the fixed maximum diameter, adn the material

%recover fixed shape characteristics
d0 = adapter.base_diam; %[m] diameter of the base
adapter.r = [d0/3, d0/2]; %[m] payload adapter r1, r2
adapter.h = d0 / 4; %[m] payload adapter height 

%recover material properties
mat = mat_switch(adapter.mat_id);

%compute adapter characteristics
[adapter.m , adapter.th] = buckling_bending(adapter, loads, mat);

end

function [fairing, loads] = fairing_fun1( m_pay, fairing, loads, plotcase)

% computes the shape and mass of the fairing.
% based on cubesats volumetric density constraints, takes as input also the
% aerodynamic loads.
% fairing is simplified as a cylinder and a cone.

%cubesats characteristics:
vol_den = 1.33 * 1e3; %[kg/m^3] volumetric density of cubesats
a = 0.1; %[m] cubesat unit edge length 
b = 0.05; %[m] margin distance between payload and fairing

%recover fixed shape characteristics
d0 = fairing.base_diam; %[m] diameter of the base
L_nose = d0 * 2; %[m] length of the conical nose

%recover loads:
nx = loads.nx; %longitudinal load factor [-]
nz = loads.nz; %transversal load factor [-]
K = loads.K; %factor of safety [-]

%compute payload volume (considering cubesats sizes)
% V = m_pay_max / vol_den; %[m^3] volume of max payload
V = m_pay / vol_den; %[m^3] volume of max payload
V_real = m_pay / vol_den; %[m^3] volume of real payload

%find usable base for payload
y = 0 : a : d0/2; 
N = length(y) - 1; 
base = 0;
x = zeros(N, 1);
for j = 1 : N
    x(j) = floor( sqrt( 0.25*d0^2 - y(j)^2 - b ) / a ) * a; %[m] at y(i) level, how much space i have to fill with cubesats
    base = base + 4 * a * x(j); %[m^2] update base value 
end
L_min = V / base; %[m] height of the ammissible payload
L_real = V_real / base; %[m] height of the real payload 

%get adapter characteristics
h_ad = d0 / 4;
r1_ad = d0 / 3;
r2_ad = d0 / 2;

%find height of the cylindrical part of the fairing
% L1 = L_min + b + h_ad; %[m] 
L1 = 0; %[m] as per EDBERG-COSTA

%surface of the fairing
S_cyl = L1 * pi * d0; %[m^2] cylinder surface
S_cone = pi * sqrt( d0^2 / 4 + L_nose^2 ) * d0/2; %[m^2] conical surface
fairing.S = S_cyl + S_cone; %[m^2] total surface

%compute mass (from "Launch and Entry Vehicle Design, Univ. Maryland, D.L.Akin")
fairing.m = 4.95 * fairing.S ^ 1.15;%9.89 * fairing.S; % %[kg]  

%retrieve fairing parameters
% fairing.M = nose.M + cyl.M; %[kg] total mass
fairing.L = L1 + L_nose; %[m] total height
fairing.d = d0; %[m] maximum diameter
% fairing.th_cone = nose.th; %[m] thickness of the nose
% fairing.th_cyl = cyl.th; %[m] thickness of the cylinder
fairing.V = V; %[m^3] maximum volume for the payload

%center of mass 
h_cyl = L1 / 2;
h_cone = (1/3) * L_nose + L1;
m_cyl = fairing.m * S_cyl / fairing.S;
m_cone = fairing.m * S_cone / fairing.S;
h_CG = ( m_cone * h_cone + m_cyl * h_cyl ) / fairing.m;
fairing.hCG = h_CG; %[m] barycenter height measured from the base of the fairing

%recover loads acting on base
[load, hCG] = bending_load_weight(loads, m_cone, m_cyl, h_cyl, L1, (1/3) * L_nose);
loads.F = load.F; %[N] shear force on the base of the fairing
fairing.h_m = hCG; %[m] barycenter height measured from the base of the fairing

%plot of the fairing
if nargin > 3
    h0 = fairing.h0;
    X = [d0/2, d0/2, 0, -d0/2, -d0/2];
    Y = h0 + [0, L1, fairing.L, L1, 0];
    plot(X, Y, 'k', DisplayName='true'); grid on; hold on; axis equal;%fairing
    X_a = [0, r2_ad, r1_ad, -r1_ad, -r2_ad, 0];
    Y_a = h0 + [0, 0, h_ad, h_ad, 0, 0];
    plot(X_a, Y_a, '--r', DisplayName='true'); %adapter
    X_p = [0, x(1), x(1), -x(1), -x(1), 0];
    Y_p = h0 + [h_ad, h_ad, h_ad + L_real, h_ad + L_real, h_ad, h_ad]; 
    plot(X_p, Y_p, 'b', DisplayName='true'); %payload
    X_r = [0, 0.6, 0.6, -0.6, -0.6, 0];
    Y_r = [0, 0, h0, h0, 0, 0];
    plot(X_r, Y_r, '--k', DisplayName='true');
    legend('Fairing', 'Adapter', 'Payload','Hypothetical Rocket', 'interpreter', 'latex');
    xlabel('x [m]', 'Interpreter','latex');
    ylabel('y [m]', 'Interpreter','latex');
end
end

function [loads, hCG] = bending_load_weight(loads, m1, m2, hCG1, h1, hCG2)
%
%               ^ F
%        ^ F1   |  ^ F2
%    ____|______|__|______|\
%    |___m2____||__m1_____|\  + M
%               °hCG      |\
%                      <--|
%                     h   |
%

%recover loads parameters
nz = loads.nz; %[-] transversal acceleration load factor
% L = loads.lift; %[N] vector of the concentrated lift forces [nx1]
% x = loads.X_lift; %[m] vector of the contentrated lift force positions [nx1]
F_m = (m1+m2) * 9.81 * nz; %[N] shear force acting on the base of the structure

%compute arm length of the resultant application
hCG = ( (h1+hCG2)*m2 + hCG1*m1 ) / (m1+m2); %[m] 

%compute the bending moment acting on the base of the structure
M = F_m * hCG ;%+ x * L'; %[Nm] resulting bending moment

%loads
loads.M = M;
loads.F = F_m ;%+ sum( L );

end

