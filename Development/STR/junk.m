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