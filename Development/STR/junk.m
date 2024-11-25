%tanks shape definition
e = sqrt( 1 - 1 / AR^2 );         %eccentricity of oblate part [-]
R_int = @(t) (diam - 2*t) / 2;    %internal radius of the tank [m]
V_obl = @(t) (4/3)*pi*R_int(t)^3 / AR; %volume of the two oblate parts [m^3]

%spherical tank hypothesis
R_sphere_lox = ( 3 * vlox / ( 4*pi ) )^(1/3); %[m] hypothetical radius of a spherical tank
R_sphere_rp1 = ( 3 * vrp1 / ( 4*pi ) )^(1/3); %[m] hypothetical radius of a spherical tank

if R_sphere_lox < 0.98 * diam/2
    y_lox = @(t) 2*R_sphere_lox;%fluid level inside the tank [m]
    l_lox = @(t) y_lox(t) + 2*t;%height of the tank [m]
    S_lox = @(t) 4 * pi * R_sphere_lox^2;%surface of tank [m^2]
else 
    V_cyl_lox = @(t) vlox - V_obl(t); %volume of the cylindrical part [m^3]
    h_cyl_lox = @(t) V_cyl_lox(t) / (pi*R_int(t)^2); %height of cylindrical part [m]
    y_lox = @(t) h_cyl_lox(t) + 2*R_int(t)/AR;%fluid level inside the tank [m]
    l_lox = @(t) y_lox(t) + 2*t;%height of the tank [m]
    S_cyl_lox = @(t) 2*pi*R_int(t)* h_cyl_lox(t); %surface of cylindrical part [m^2]
    S_obl_lox = @(t) 2*pi * R_int(t)^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
    S_lox = @(t) S_obl_lox(t) + S_cyl_lox(t); %surface of the lox tank [m^2]
end

if R_sphere_rp1 < 0.98 * diam/2
    y_rp1 = @(t) 2*R_sphere_rp1;%fluid level inside the tank [m]
    l_rp1 = @(t) y_rp1(t) + 2*t;%height of the tank [m]
    S_rp1 = @(t) 4 * pi * R_sphere_rp1^2;%surface of tank [m^2]
else
    V_cyl_rp1 = @(t) vrp1 - V_obl(t); %volume of the cylindrical part [m^3]
    h_cyl_rp1 = @(t) V_cyl_rp1(t) / (pi*R_int(t)^2); %height of cylindrical part [m]
    y_rp1 = @(t) h_cyl_rp1(t) + 2*R_int(t)/AR;%fluid level inside the tank [m]
    l_rp1 = @(t) y_rp1(t) + 2*t;%height of the tank [m]
    S_cyl_rp1 = @(t) 2*pi*R_int(t)* h_cyl_rp1(t); %surface of cylindrical part [m^2]
    S_obl_rp1 = @(t) 2*pi * R_int(t)^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
    S_rp1 = @(t) S_obl_rp1(t) + S_cyl_rp1(t); %surface of the rp1 tank [m^2]
end


%pressure at base with longitudinal acceleration (Stevin's law)
p_lox = @(t) p + y_lox(t) * rholox * acc; %[Pa] pressure at bottom of tank during acceleration
p_rp1 = @(t) p + y_rp1(t) * rhorp1 * acc; %[Pa] pressure at bottom of tank during acceleration

%thickness of tanks
f_lox = @(t) t - (diam - 2*t)*p_lox(t)/( 2*sy );
f_rp1 = @(t) t - (diam - 2*t)*p_rp1(t)/( 2*sy );
t_lox = fzero(f_lox, [t_min, 1]); %[m] lox tank thickness
t_rp1 = fzero(f_rp1, [t_min, 1]); %[m] rp1 tank thickness
t.lox = t_lox;
t.rp1 = t_rp1;

%check on manufacturability
if t_lox < t_min
    t_lox = t_min;
elseif t_rp1 < t_min
    t_rp1 = t_min;
end

%height of tanks
h.tank_lox = l_lox(t_lox); %[m] height of tank
h.tank_rp1 = l_rp1(t_rp1); %[m] height of tank
h.tot = h.tank_lox + h.tank_rp1; %[m] total height of tanks together

%surfaces estimation
Slox = S_lox(t_lox);%surface of the lox tank [m^2]
Srp1 = S_rp1(t_rp1);%surface of the rp1 tank [m^2]

%%
clc;
clear;
close all;
% Initial constants
m_pay = 5000;            % Payload mass in kg
dv = 3000;               % Delta-v in m/s (for each stage)
OF = 2.5;                % Oxidizer to Fuel ratio
Is = 300;                % Specific impulse in seconds for both stages
eps_start = 0.05;        % Initial guess for the structural mass index (eps)
m_motors1 = 1000;        % Motor mass for stage 1 in kg
m_motors2 = 800;         % Motor mass for stage 2 in kg
d = 2.5;                 % Diameter in meters
AR = 6;                  % Aspect ratio for the stages
loads = [1.2, 1.5];      % Load factors for stages
mat1 = 'Aluminum';       % Material for stage 1
mat2 = 'Titanium';       % Material for stage 2
press1 = 10;             % Pressure for stage 1 in MPa
press2 = 15;             % Pressure for stage 2 in MPa

% Tolerance and iteration parameters
tol = 1e-4;              % Desired tolerance for eps difference
max_iter = 1000;         % Maximum number of iterations
iter = 0;                % Initialize iteration counter
learning_rate = 1e-3;    % Step size for updating eps

% Initialize other variables
eps_end = [0.05; 0.05];  % Initial eps_end guess (can be set based on system)
eps_diff = abs(eps_start - eps_end); % Initial difference

% Initialize GLOM struct (you can define this struct as per your system)
GLOM.m_stag = zeros(2,1); % Placeholder
GLOM.m_tot = zeros(2,1);  % Placeholder
GLOM.m_prop = zeros(2,1); % Placeholder

% Initialize M struct for stage mass indexes (again, set based on system)
M.eps1 = eps_start;
M.eps2 = eps_start;

% Begin iterative process
while max(eps_diff) > tol && iter < max_iter
    iter = iter + 1; % Increment iteration count
    
    % Call TANDEM function to calculate masses and propellants
    [GLOM.m_stag, GLOM.m_tot, GLOM.m_prop] = TANDEM(Is, eps_start, dv, m_pay, 0);
    
    % Recalculate the first stage mass properties (M1)
    M1.rp1 = GLOM.m_prop(1) * 1 / (1 + OF);
    M1.lox = GLOM.m_prop(1) * OF / (1 + OF);
    M1.prop = M1.lox + M1.rp1;
    M1.motor = m_motors1;
    
    % Calculate inert mass for stage 1
    [M1, h1, th1] = inert_mass(M1, d, AR, loads, mat1, press1);
    
    % Recalculate the second stage mass properties (M2)
    M2.rp1 = GLOM.m_prop(2) * 1 / (1 + OF);
    M2.lox = GLOM.m_prop(2) * OF / (1 + OF);
    M2.prop = M2.lox + M2.rp1;
    M2.motor = m_motors2;
    
    % Calculate inert mass for stage 2
    [M2, h2, th2] = inert_mass(M2, d, AR, loads, mat2, press2);
    
    % Update the eps1 and eps2 values in M struct based on new mass calculations
    M.eps1 = M1.eps;
    M.eps2 = M2.eps;
    
    % Calculate new eps_end based on updated stage properties
    eps_end = [M.eps1; M.eps2];
    
    % Calculate the difference between the current eps and eps_end
    eps_diff = abs(eps_start - eps_end);
    
    % Update eps using a simple gradient descent-like step
    if eps_end(1) > eps_start
        eps_start = eps_start + learning_rate;  % Increase eps if eps_end is greater
    else
        eps_start = eps_start - learning_rate;  % Decrease eps if eps_end is smaller
    end
    
    % Print debug information for each iteration (optional)
    disp(['Iteration ', num2str(iter), ' - eps_diff = ', num2str(max(eps_diff))]);
end

% Check if the solution converged
if max(eps_diff) <= tol
    disp(['Convergence achieved in ', num2str(iter), ' iterations.']);
else
    disp('Maximum iterations reached without convergence.');
end

%%


clc;
clear;
close all;

%data from other departments:
Is = [311; 311]; %[s] stages Is
dv = 8.5; %[km/s] required dv
M.pay = 250; %[kg] payload mass
OF = 2.58; %[-] Ox/Fu ratio for LOX-RP1
n = 5; %[-] load factor of longitudinal acceleration
t = 5; %[-] load factor of transversal acceleration
diam = 1.2; %[m] external diameter
AR = sqrt(3);%sqrt(3); %aspect ratio of oblate domes [-]
loads.n = n; %longitudinal acceleration [-]
loads.t = t;%transversal load factor [-]
loads.K = 1.5; %loads resistance safety factor [-]

%stage 1
M1.OF = OF;%[-] Ox/Fu ratio
M1.motor = 315; %[kg] only motor, pumps and batteries (electron - rutherford motor) %pump-fed
M1.fairing = 100; %[kg] fairing of the first stage is nonexistent
M1.rhorp1 = 807;  %[kg/m^3] density of rp1
M1.rholox = 1140; %[kg/m^3] density of lox
mat1 = 5; % 1 for Ti, 2 for Al 2XXX, 3 for Steel, 4 for Carbon Fiber, 5 for Al 7XXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% in future versions can be optimized the material selection in function
press1 = 2; % 0 for unpressurized, 1 for pressure-fed, 2 for pump-fed, 3 for blowdown

%stage 2
M2.OF = OF;%[-] Ox/Fu ratio
M2.motor = 45; %[kg] only motor, pumps and batteries (electron - rutherford motor) %pump-fed
M2.fairing = 30; %[kg] fairing of the second stage (31.8 / 31.9)
M2.rhorp1 = 807;  %[kg/m^3] density of rp1
M2.rholox = 1140; %[kg/m^3] density of lox
mat2 = 1; % 1 for Ti, 2 for Al 2XXX, 3 for Steel, 4 for Carbon Fiber, 5 for Al 7XXX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% in future versions can be optimized the material selection in function
press2 = 2; % 0 for unpressurized, 1 for pressure-fed, 2 for pump-fed, 3 for blowdown

%first guesses:
eps0 = [0.1; 0.2]; %[0.06; 0.2]; %[-] stages structural mass indexes



shape.r = [diam/3, diam/2];%[m]
shape.h = 0.4; %[m]
loads.m = 250; %[kg]

mat.rho = 4500; %[kg/m^3]
mat.t_min = 0.5 * 1e-3; %[m] minimum thickness for manufacturability
mat.E = 110 * 1e9; %[Pa] young modulus
mat.sy = 900 * 1e6; %[Pa] tensile yield stress
mat.su = 950 * 1e6; %[Pa] tensile ultimate stress
mat.nu = 0.34; %[-] Poisson's ratio

[m, th] = buckling(shape, loads, mat, press1);


%%

    % if R_sphere_lox < R_int %that is, if both tanks are spherical
    %     shape2.r = [R_sphere_rp1, R_sphere_lox, R_int];
    % else %that is, if rp1 is spherical but lox is cylindrical
    %     shape2.r = [R_sphere_rp1, R_sphere_lox, R_int]; %for simplicity we take the same dimensions of the "both-spherical" case
    % end


%%

%COMPUTE MASSES AND THICKNESSES:

%conical nose:
nose.r = [0, d0/2]; %nose radiuses [m]
nose.h = L_nose; %nose lenght [m]
nose_loads.m = 0; %sustained mass [kg]
nose_loads.n = n; %longitudinal load factor [-]
nose_loads.K = K; %factor of safety [-]
nose_loads.F_drag = F_aero; %aerodynamic drag force [N]
[nose.M, nose.th] = buckling(nose, nose_loads, mat, 0);

%cylindrical body:
cyl.r = d0/2; %cylinder radius [m]
cyl.h = L1; %cylinder height [m]
cyl_loads.m = nose.M; %sustained mass [kg]
cyl_loads.n = n; %longitudinal load factor [-]
cyl_loads.K = K; %factor of safety [-]
cyl_loads.F_drag = F_aero; %aerodynamic drag force [N]
[cyl.M, cyl.th] = buckling(cyl, cyl_loads, mat, 0);


%%

if id == 4 %composites need different formulas

    if p == 0 %FOR CONNECTORS ( p = 0 )
        if length(r) == 1 %cylindrical connector
           
        else %CONICAL CONNECTOR
         
        end
    else %FOR CYLINDRICAL TANKS ( p =\= 0, size(r)==1 )
    
    end
else
    if p == 0 %FOR CONNECTORS ( p = 0 )
        if length(r) == 1 %CYLINDRICAL CONNECTOR
            %compute knockdown factor:
            phi = @(th) (1/16) * sqrt(r/th);
            gP = @(th) 1 - 0.901 * ( 1 - exp( -phi(th) ) );
            gM = @(th) 1 - 0.731 * ( 1 - exp( -phi(th) ) );

            %get the wall flexural stiffness per unit width:
            D = @(th) E * th^3 / ( 12 * (1-nu^2) );

            %critical loads (axial compression and bending):
            Pcr = @(th) k_x(nu, L, r, th, gP) * 2 * pi^3 * D(th) * r / L^2;
            Mcr = @(th) k_x(nu, L, r, th, gM)   *   pi^3 * D(th) * r^2 / L^2;
            %simplified expression for Mcr : Mcr = @(th) 0.987 * E * th^2 / sqrt(1-nu^2);

            %minimum thickness to sustain only Pcr:
            t_l = fzero( K * F_load - Pcr , [r/1500, 1]);

        else %CONICAL CONNECTOR 
    
        end
    else %FOR CYLINDRICAL TANKS ( p =\= 0, size(r)==1 )
        %compute knockdown factor:
        phi = @(th) (1/16) * sqrt(r/th);
        gP = @(th) 1 - 0.901 * ( 1 - exp( -phi(th) ) );
        gM = @(th) 1 - 0.731 * ( 1 - exp( -phi(th) ) );
        dg = @(th) d_gamma(p, E, r, th); 

        %critical loads (axial compression and bending):
        Pcr = @(th) 2*pi    * E * th^2 * ( gP(th) / sqrt( 3 * (1-nu^2) ) + dg(th) ) + p * pi * r(1)^2; 
        Mcr = @(th) pi*r(1) * E * th^2 * ( gM(th) / sqrt( 3 * (1-nu^2) ) + dg(th) ) + p * pi * r(1)^2 * k1;

        %minimum thickness to sustain only Pcr:
        t_l = fzero( K * F_load - Pcr , [r/1500, 1]);
    end
end




%%


% f(eps)
% f(1)
% th = fzero( f , t_min2); %[m]
% % x = 0:1e-5:0.01;
% % fx = abs(f(x));
% % [~, min_pos] = min(fx);
% % th = x(min_pos);%[m]
% %check on manufacturability (th>t_min) and on validity of the NASA model (th>th_lim = r/1500):
% %th = max(th, t_min);
% th = max([th, t_min1, t_min2]); %[m] 

%%

% %adapter 
% adapter.m = 0.0755 * M.pay + 50; %[kg] estimated mass from Edberg-Costa

% %adapter 
% loads_a = loads;
% loads_a.m = M.pay;%sustained mass [kg]
% loads_a.F_drag = 0; %aerodynamic load is null for the payload adapter [N]
% loads_a.p = 0; %internal pressure [Pa]
% [adapter] = adapter_fun1(adapter, loads_a);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get adapter characteristics
adapter_loads = loads;
adapter_loads.m = m_pay; %sustained mass [kg]
adapter_loads.F_drag = 0; %aerodynamic drag force [N] is null for the adapter
adapter_loads.p = 0; %internal pressure [Pa]
adapter.mat_id = fairing.mat_id; %in reality, we just need the height of the adapter for this step, therefore we can put any material, since the shape is fixed
adapter.base_diam = fairing.base_diam; %[m] diameter of the base
adapter = adapter_fun1(adapter, adapter_loads);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

