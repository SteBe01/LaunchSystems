function [Tank] = tank_mass(M, diam, AR, loads, mat, pressure_type)

% considers thickness equal along the whole tank.
% safety factor to be defined.
% the volume is the volume of propellant to be contained: you cannot use
% this function to evaluate blowdown architectures.

%propellant masses
OF = M.OF;%[-] Ox/Fu ratio
Tank.OF = OF;
mlox = M.prop * OF / (1+OF);%[kg] mass of lox
mrp1 = M.prop * 1  / (1+OF);%[kg] mass of rp1
Tank.OX.M_lox = mlox;
Tank.FU.M_fuel = mrp1;
% mlox = M.lox; %[kg]
% mrp1 = M.rp1; %[kg]
% M.prop = mrp1 + mlox; %[kg] total propellant mass

%propellant densities
rholox = M.rholox; %[kg/m^3]
rhorp1 = M.rhorp1; %[kg/m^3]
Tank.OX.rho_lox = rholox;
Tank.FU.rho_fuel = rhorp1;
%propellant volumes
vlox = mlox / rholox; %[m^3]
vrp1 = mrp1 / rhorp1; %[m^3]
Tank.OX.V_lox = vlox;
Tank.FU.V_fuel = vrp1;

%acceleration
acc = loads.acc;

switch pressure_type 
    case 0 % case for unpressurized vessel
        MEOP = 0;
    case 1 %case for pressure-fed
        MEOP = 2.8 * 1e6; %[Pa] internal tank pressure
    case 2 %case for pump-fed
        MEOP = 300 * 1e3; %[Pa] internal tank pressure
    case 3 % case for blowdown
        MEOP = 50 * 1e6; %[Pa] internal tank pressure
end

switch mat 
    case 1 % Ti6Al4V
        rho = 4500; %[kg/m^3]
        t_min = 0.5 * 1e-3; %[m] minimum thickness for manufacturability
        E = 110 * 1e9; %[Pa] young modulus
        sy = 900 * 1e6; %[Pa] tensile yield stress
        su = 950 * 1e6; %[Pa] tensile ultimate stress
    case 2 % Al 2XXX
        rho = 2700; %[kg/m^3]
        t_min = 0.5 * 1e-3; %[m] minimum thickness for manufacturability
        E = 70 * 1e9; %[Pa] young modulus
        sy = 290 * 1e6; %[Pa] tensile yield stress
        su = 390 * 1e6; %[Pa] tensile ultimate stress
    case 3 % Steel
        rho = 7800; %[kg/m^3]
        t_min = 0.25 * 1e-3; %[m] minimum thickness for manufacturability
        E = 200 * 1e9; %[Pa] young modulus
        sy = 350 * 1e6; %[Pa] tensile yield stress
        su = 420 * 1e6; %[Pa] tensile ultimate stress
    case 4 % Carbon fiber
        rho = 1800; %[kg/m^3]
        t_min = 0.90 * 1e-3; %[m] minimum thickness for manufacturability
        E = 250 * 1e9; %[Pa] young modulus
        sy = 350 * 1e6; %[Pa] tensile yield stress
        su = sy; %[Pa] tensile ultimate stress
end

Tank.rho = rho;
Tank.t_min = t_min;
Tank.E = E;
Tank.sy = sy;
Tank.su = su;

%correction factor
Tank.MEOP = MEOP;
Km = 1.1; 
Tank.Km = Km;
Kp = 1.5;
Tank.Kp = Kp;
MDP = MEOP * Km * Kp;
jproof = 1.25;
jburst = 1.5;
p = MDP * jburst;
Tank.jburst = jburst; % Burst pressure is used!!!
%Tank.jproof = jproof;
Tank.p = p;
%p = MEOP;


%tanks shape definition
e = sqrt( 1 - 1 / AR^2 );         %eccentricity of oblate part [-]
R_int = diam / 2;    %internal radius of the tank [m]
V_obl = (4/3)*pi*R_int^3 / AR; %volume of the two oblate parts [m^3]
%spherical tank hypothesis
R_sphere_lox = ( 3 * vlox / ( 4*pi ) )^(1/3); %[m] hypothetical radius of a spherical tank
R_sphere_rp1 = ( 3 * vrp1 / ( 4*pi ) )^(1/3); %[m] hypothetical radius of a spherical tank

if R_sphere_lox < 0.98 * diam/2
    y_lox = 2*R_sphere_lox;%fluid level inside the tank [m] (hyp: full tank)
    l_lox = @(t) y_lox + 2*t;%height of the tank [m] (sum the two thicknesses)
    S_lox = 4 * pi * R_sphere_lox^2;%surface of tank [m^2]
    Tank.OX.y_lox = y_lox;
    Tank.OX.S_lox = S_lox;
    Tank.OX.R_sphere_lox = R_sphere_lox;
    Tank.OX.FLAG = 'Sphere';
else 
    V_cyl_lox = vlox - V_obl; %volume of the cylindrical part [m^3]
    h_cyl_lox = V_cyl_lox / (pi*R_int^2); %height of cylindrical part [m]
    y_lox = h_cyl_lox + 2*R_int/AR;%fluid level inside the tank [m]
    l_lox = @(t) y_lox + 2*t;%height of the tank [m]
    S_cyl_lox = 2*pi*R_int* h_cyl_lox; %surface of cylindrical part [m^2]
    S_obl_lox = 2*pi * R_int^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
    S_lox = S_obl_lox + S_cyl_lox; %surface of the lox tank [m^2]
    Tank.OX.V_oblate = V_obl;
    Tank.OX.S_lox = S_lox;
    Tank.OX.y_lox = y_lox;
    Tank.OX.R_cyl_lox = R_int;
    Tank.OX.FLAG = 'Cyl';
    Tank.OX.H_cyl = h_cyl_lox;
    Tank.OX.H_dome = R_int/AR;
end

if R_sphere_rp1 < 0.98 * diam/2
    y_rp1 = 2*R_sphere_rp1;%fluid level inside the tank [m]
    l_rp1 = @(t) y_rp1 + 2*t;%height of the tank [m]
    S_rp1 = 4 * pi * R_sphere_rp1^2;%surface of tank [m^2]
     Tank.FU.y_fuel= y_rp1;
    Tank.FU.S_fuel = S_rp1;
    Tank.FU.R_sphere_fuel = R_sphere_rp1;
    Tank.FU.FLAG = 'Sphere';
else
    V_cyl_rp1 = vrp1 - V_obl; %volume of the cylindrical part [m^3]
    h_cyl_rp1 = V_cyl_rp1 / (pi*R_int^2); %height of cylindrical part [m]
    y_rp1 = h_cyl_rp1 + 2*R_int/AR;%fluid level inside the tank [m]
    l_rp1 = @(t) y_rp1 + 2*t;%height of the tank [m]
    S_cyl_rp1 = 2*pi*R_int* h_cyl_rp1; %surface of cylindrical part [m^2]
    S_obl_rp1 = 2*pi * R_int^2 * ( 1 + (1/(e*AR^2)) * atanh( e ) ); %surface of the oblate parts together [m^2]
    S_rp1 = S_obl_rp1 + S_cyl_rp1; %surface of the rp1 tank [m^2]
    Tank.FU.S_fuel = S_rp1;
    Tank.FU.y_fuel = y_rp1;
    Tank.FU.R_cyl_fuel = R_int;
    Tank.FU.FLAG = 'Cyl';
    Tank.FU.H_cyl = h_cyl_rp1;
    Tank.FU.H_dome = R_int/AR;
end


%pressure at base with longitudinal acceleration (Stevin's law)
p_lox = p + y_lox * rholox * acc; %[Pa] pressure at bottom of tank during acceleration
p_rp1 = p + y_rp1 * rhorp1 * acc; %[Pa] pressure at bottom of tank during acceleration
Tank.p_lox_bottom = p_lox;
Tank.p_fuel_bottom = p_rp1;
%thickness of tanks
t_lox = diam*p_lox/( 2*sy ); %[m] lox tank thickness
t_rp1 = diam*p_rp1/( 2*sy ); %[m] rp1 tank thickness
t.lox = t_lox;
t.rp1 = t_rp1;

%check on manufacturability
if t_lox < t_min
    t_lox = t_min;
elseif t_rp1 < t_min
    t_rp1 = t_min;
end

Tank.OX.t_lox = t_lox;
Tank.FU.t_fuel = t_rp1;

%height of tanks
h.tank_lox = l_lox(t_lox); %[m] height of tank
h.tank_rp1 = l_rp1(t_rp1); %[m] height of tank
h.tot = h.tank_lox + h.tank_rp1; %[m] total height of tanks together 
%                                     (one on top of the other)

Tank.OX.tank_lox = h.tank_lox;
Tank.FU.tank_fuel = h.tank_rp1;
Tank.h_tot = h.tot;

if R_sphere_lox < 0.98 * diam/2
 
m_sphere_lox = 2*pi*((R_sphere_lox)^3) * (p_lox/(sy/rho));
M_tot_tank_lox = m_sphere_lox;
Tank.OX.M_tot_tank_lox = M_tot_tank_lox;
else 
    
M_spherical_cap_lox = rho * ((pi*((2*R_int)^3))/(4*sy)) * (p + (0.5*rholox*acc*((h.tank_lox/(2*R_int)) - (1/6))));
Tank.OX.M_spherical_cap_lox = M_spherical_cap_lox;
M_cyl_lox = rho * ((pi*((2*R_int)^3))/(4*sy)) * ((2*p*((h.tank_lox/(2*R_int))-1)) + (rholox*acc*2*R_int*((((h.tank_lox/(2*R_int))-1))^2)) );
Tank.OX.M_cyl_lox = M_cyl_lox;
M_tot_tank_lox = (2*M_spherical_cap_lox) + M_cyl_lox;
Tank.OX.M_tot_tank_lox = M_tot_tank_lox;
end

if R_sphere_rp1 < 0.98 * diam/2
 
m_sphere_rp1= 2*pi*((R_sphere_rp1)^3) * (p_rp1/(sy/rho));
M_tot_tank_fuel = m_sphere_rp1;
Tank.FU.M_tot_tank_fuel = M_tot_tank_fuel;
else 
    M_spherical_cap_fuel = rho * ((pi*((2*R_int)^3))/(4*sy)) * (p + (0.5*rhorp1*acc*((h.tank_rp1/(2*R_int)) - (1/6))));
    Tank.FU.M_spherical_cap_fuel = M_spherical_cap_fuel;
M_cyl_fuel = rho * ((pi*((2*R_int)^3))/(4*sy)) * ((2*p*((h.tank_rp1/(2*R_int))-1)) + (rhorp1*acc*2*R_int*((((h.tank_rp1/(2*R_int))-1))^2)) );
Tank.FU.M_cyl_fuel = M_cyl_fuel;
M_tot_tank_fuel = (2*M_spherical_cap_fuel) + M_cyl_fuel;
Tank.FU.M_tot_tank_fuel = M_tot_tank_fuel;
end

Tank.M_tot_struct = M_tot_tank_lox + M_tot_tank_fuel;

Tank.M_total = Tank.M_tot_struct + Tank.OX.M_lox + Tank.FU.M_fuel;

Tank.M_prop_tot =  Tank.OX.M_lox + Tank.FU.M_fuel;

end