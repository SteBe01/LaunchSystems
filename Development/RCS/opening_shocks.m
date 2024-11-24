%% Parachute design
clc; clear
m = 1090;              % [kg]
Cd0 = 0.8;              % [-] from Paper on ringsail (modified)
h_reef = 10;        % [km]
h_disreef = 3;      % [km]
v_end = 7;             %[m/s] design touchdown velocity
rho0 = 1.225;
g0 = 9.81;

rho_reef = density_model(h_reef);    % [kg/m^3]
rho_disreef = density_model(h_disreef);    % [kg/m^3]
g_reef = gravity(h_reef);
g_disreef = gravity(h_disreef);

S_disreef = (m*g0)./(0.5*rho0*v_end^2*Cd); % surface of main parachute (disreefed)

% Shock Main parachute

perc = 0.15; %assuming 30 percent of initial surface
S_reefed = perc*S_disreef;

v_reef = sqrt((m*g_reef)/(0.5*rho_reef*Cd*S_reefed));
v_drogue = 80;   % assuming 

A_reef = ballistic_par(S_disreef,rho_reef,v_drogue,Cd,m,perc,1);
% A_reef = 0.6340--> X1 = 0.35
X1_reef = 0.35; % opening force reduction factor
F_reef = Cd*S_reefed*0.5*rho_reef*v_drogue^2* X1_reef


A_disreef = ballistic_par(S_disreef,rho_disreef,v_reef,Cd,m,perc,2);
% A_disreefed = 0.1896 --> X1 = 0.2
X1_disreef = 0.2;
F_disreef = Cd*S_disreef*0.5*rho_reef*v_reef^2*X1_disreef



% define pilot and drogue chutes area basing on literature percentage of
% main parachute (Knacke)

%%
S_pilot = 0.005* S_end;
S_drogue = 0.3*S_end;
v_drogue = 80;
g_drogue = gravity(9);
Cd_drogue = (m*g_drogue)./(0.5*rho*v_drogue^2*S_drogue);

% opening shocks of main parachute
%v_drogue = sqrt( )


% Pilot: Mach ? --> mach 1.5

% Drogue: Mach 1.5 --> mach 0.3

function rho = density_model(h)

%input  h: altitude in [km]

h0_vect = [0 25 30 40 50 60 70 80 90 100 110 120 130 140 150 ...
        180 200 250 300 350 400 450 500 600 700 800 900 1000]';
    rho0_vect = [1.225 3.899*1e-2 1.774*1e-2 3.972*1e-3 1.057*1e-3 ...
        3.206*1e-4 8.770*1e-5 1.905*1e-5 3.396*1e-6 5.297*1e-7 9.661*1e-8 ...
        2.438*1e-8 8.484*1e-9 3.845*1e-9 2.070*1e-9 5.464*1e-10 ...
        2.789*1e-10 7.248*1e-11 2.418*1e-11 9.158*1e-12 3.725*1e-12 ...
        1.585*1e-12 6.967*1e-13 1.454*1e-13 3.614*1e-14 1.170*1e-14 ...
        5.245*1e-15 3.019*1e-15]';
    H_vect = [7.249 6.349 6.682 7.554 8.382 7.714 6.549 5.799 5.382 5.877 ...
        7.263 9.473 12.636 16.149 22.523 29.740 37.105 45.546 53.628 ...
        53.298 58.515 60.828 63.822 71.835 88.667 124.64 181.05 268.00]';

if h < 0
    rho0 =0;
    h0 = 0;                            % [km]
    H = 0;                              % [km]
else
    idx = sum(h0_vect < h);

    rho0 = rho0_vect(idx);
    h0 = h0_vect(idx);                            % [km]
    H = H_vect(idx);                              % [km]
end
    rho = rho0 * exp(-(h-h0)/H);                % [kg/m^3]

end

function g = gravity(h)
% h = altitude [km]

Re = 6371; % [km]
r = (h+Re)*10^3;  % [m]
G = 6.6725985e-11;  % [N m2 / kg] 
M = 5.972e24;  % kg
g = G*M/r^2;

end


function A = ballistic_par(S0,rho,vs,Cd,M,perc,reefing)
% S surfacee [m^2]
% rho density at opening altitude [kg/m^3]
% v velocity of the body at opening [m/s]
% Cd 
% m mass [kg]
% perc  percentage of reefed area



D0 = sqrt(4*S0/pi);  %[ft]

%rho = rho* 2.20462 /(3.28084)^3;            %[lb/ft^3]


if reefing == 1             %reefing case
    n = 8;    % 7-8
    tf = n*D0/vs*sqrt(perc);       % filling time (canopy inflation time)
    A = 2*M/(Cd*S0*perc*rho*vs*tf);
elseif reefing == 2         % disreefing case
     n = 2;    
     tf = n*D0/vs*sqrt(1-perc);       % filling time (canopy inflation time)
     A = 2*M/(Cd*S0*rho*vs*tf);
end

end
