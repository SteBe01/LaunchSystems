



kx = @(th) k_x(0.3, 1, 0.6, th, 0.41);

kx(0.0001)
%% Functions

function [th, XY] = buckling_bending(shape, load, mat, plotcase)

% based on NASA papers in shared folder (SP-8007 & SP-8019)
% computes connectors masses, heights and thicknesses to sustain
% compression loads and avoid buckling effect

%constants:
g = 9.81; %[m/s^2] gravitational acceleration

%recover loads:
m = load.m; %sustained mass [kg]
n = load.n; %longitudinal load factor [-]
K = load.K; %factor of safety [-]
F_aero = load.F_drag; %aerodynamic drag force [N]
p = load.p; %internal pressure [Pa]
M_exp = load.M_exp; %expected flessional 

%recover material characteristics:
id = mat.ID; %[-] ID of the material: 1 for Ti, 2 for Al 2XXX, 3 for Steel, 4 for Carbon Fiber, 5 for Al 7XXX
E = mat.E; %[Pa] Young modulus
rho = mat.rho; %[kg/m^3]
t_min = mat.t_min; %[m] minimum thickness for manufacturability 
sy = mat.sy; %[Pa] tensile yield stress
su = mat.su; %[Pa] tensile ultimate stress
nu = mat.nu; %[-] Poisson's ratio

%recover dimensions:
r = shape.r;
h = shape.h; %distance between the base of the domes
if length(r) > 1 %trucated-cone
    %recover cone shape characteristics
    alpha = asin( ( r(2) - r(1) ) / h );
    L = sqrt( h^2 - (r(2) - r(1))^2 );
    l = cos(alpha) * L; %height of the shell
    rho1 = r(1);
    rho2 = r(2);
    r(2) = cos(alpha) * r(2); 
    r(1) = cos(alpha) * r(1);

    %plot
    if nargin > 4
        h0 = shape.h0; %height at which the connector is placed 
        %for the plotting
        y = h0 + [rho2*sin(alpha/2), rho2*sin(alpha/2), rho2*sin(alpha/2)+l, rho2*sin(alpha/2)+l, rho2*sin(alpha/2), rho2*sin(alpha/2)];
        XY = [0, rho2, rho1, -rho1, -rho2, 0; y];
    end
else
    %recover cylinder shape characteristics
    alpha = 0;
    L = h;
    if nargin > 4
        h0 = shape.h0; %height at which the connector is placed 
        %for the plotting
        y = h0 + [0, 0, L, L, 0, 0];
        XY = [0, r, r, -r, -r, 0; y];
    end
end

%compute sustained load:
F_load = m * n * g + F_aero; %load [N]

%get the wall flexural stiffness per unit width:
D = @(th) E * th^3 / ( 12 * (1-nu^2) );

%compute delta_gamma for pressure-increased performances:
dg = @(th) d_gamma(p, E, r, th, alpha); 

%6 SITUATIONS: p OR NON p , CONICAL/CYLINDRICAL (BUT THERE AREN'T CONICAL TANKS), ISOTROPIC/ORTHOTROPIC
switch id
    case 4 %(for CF, orthotropic expressions)

    otherwise
        if alpha == 0 %(cylindrical shape, NASA SP-8007)
            k1 = 0.8;

            %compute knockdown factor:
            phi = @(th) (1/16) * sqrt(r(1)/th);
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
Pcr = @(th) k2 * k_x(nu, L, r, th, gP) * 2 * pi^3 * D(th) * r   / L^2 +...   %if cyl and p=0 (only metals)
         + (1-k2) * ( 2*pi    * E * th^2 * ( gP(th) / sqrt( 3 * (1-nu^2) ) + dg(th) ) + p * pi * r(1)^2 ); %in any other case (only metals)
Mcr = @(th) k2 * k_x(nu, L, r, th, gM)   *   pi^3 * D(th) * r^2 / L^2 +...   %if cyl and p=0 (only metals)
         + (1-k2) * pi*r(1) * E * th^2 * ( gM(th) / sqrt( 3 * (1-nu^2) ) + dg(th) ) + p * pi * r(1)^2 * k1;%in any other case (only metals)

%relation to be satisfied in combined stress condition
f = @(th) K*F_load/Pcr(th) + K*M_exp/Mcr(th) - 1; %this must be <0 to preserve the structure

%minimum thickness to sustain only Pcr:
th = fzero( f , [r/1500, 1]); %[m]

%check on manufacturability:
th = max(th, t_min); %[m]
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

gZ = @(th) g * L^2 / (r*th) * sqrt(1-nu^2);

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



% %%
% %two cases: cylindrical of trucated-cone shells:
% if length(r) > 1 %trucated-cone
% 
%     %recover cone shape characteristics
%     alpha = asin( ( r(2) - r(1) ) / h );
%     L = sqrt( h^2 - (r(2) - r(1))^2 );
%     l = cos(alpha) * L; %height of the shell
%     r2 = cos(alpha) * r(2);
%     r1 = cos(alpha) * r(1);
% 
%     if F < 0 %load is an axial tension load
%         th = t_min;
%     else %load is an axial compression load
%         %compute the critical thickness for the loaded shell
%         gamma = 0.33;
%         th = sqrt( ( K*F*sqrt( 3*(1-nu^2) ) ) / ( 2*pi * gamma * E * cos(alpha)^2 ) ); %see page 13 of NASA paper in launch systems shared folder
%     end
%     %compute the mass
%     S = pi * L * ( r2 + r1 ); %surface of the truncated cone
%     M = S * th * rho;
% 
%     if nargin > 4
%         h0 = shape.h0; %height at which the connector is placed 
%         %for the plotting
%         y = h0 + [r(2)*sin(alpha/2), r(2)*sin(alpha/2), r(2)*sin(alpha/2)+l, r(2)*sin(alpha/2)+l, r(2)*sin(alpha/2), r(2)*sin(alpha/2)];
%         XY = [0, r(2), r(1), -r(1), -r(2), 0; y];
%     end
% else %cylindrical shell
% 
%     %recover cylinder shape characteristics
%     l = h;
% 
%     if F < 0 %load is an axial tension load
%         th = t_min;
%     else %load is an axial compression load
%         %compute the critical thickness for the loaded shell
%         s_cr = @(t) E * ( 9 * (t/r)^1.6 + 0.16 * (t/l)^1.3 ); %[Pa] critical stress
%         F_crit = @(t) 2 * pi * r * s_cr(t) * t; %[N] critical load
%         f = @(t) F_crit(t) - K * F;
%         th = fzero(f, [0, 1]);
%     end
%     %compute the mass
%     M = 2*pi * r * l * th * rho; %mass of the connector
% 
%     if nargin > 4
%         h0 = shape.h0; %height at which the connector is placed 
%         %for the plotting
%         y = h0 + [0, 0, h, h, 0, 0];
%         XY = [0, r, r, -r, -r, 0; y];
%     end
% end
% 
% end

