function [CL, CD] = file_funzione1(Mach_v, alpha_v, h)

%
% La funzione calcola Cd e Cl prendendo in input: alpha, Mach, quota, e la geometria sotto, 
% superficie dell'area portante.
% E' i lementato Van Driest II per il calcolo della resistenza a frizione del corpo
% l'ala è modellata come 4 alette (90deg una dall'altra) immerse nel flusso,

%% INPUTS:
geometry_funzione1


%% Calculation AXIAL COEFF.:

Ca = [];
Ca_w = [];
Ca_f = [];
Ca_b = [];


for i = 1:length(Mach_v)
    %Mach = Mach_v(i);

    for j = 1:length(alpha_v)

        %alpha = alpha_v(j);       

        [Ca(i,j), Ca_w(i,j), Ca_f(i,j), Ca_b(i,j)] = drag_estimation( a, b, x, ln, dn, nose_type, Mach_v(i), alpha_v(j), h, A_w);

    end
end




%% Calculation NORMAL COEFF.:

Cn = [];

for i = 1:length(Mach_v)
    %Mach = Mach_v(i);

    for j = 1:length(alpha_v)

        %alpha = alpha_v(j);       

        [Cn(i,j),Cn_Cn0_sb(i,j),Cn_Cn0_Newt(i,j), Cdn(i,j),Cn_s(i,j)] =  CN (max(a), max(b), x, ln, dn, A_b, A_r, A_p, Mach_v(i), alpha_v(j), h, phi, data,A_w);

    end
end



%% CL and CD conversion:

[CL, CD] = body2wind(Cn, Ca, alpha_v);



%% PLOTS:


% figure
% subplot(2, 2, 1)
% plot(Mach_v, Ca)
% xlabel('Mach')
% ylabel('Ca')
% legend('alpha: ',  int2str(alpha_v))
% 
% subplot(2, 2, 2)
% plot(Mach_v, Ca_w)
% xlabel('Mach')
% ylabel('Ca_w')
% 
% 
% subplot(2, 2, 3)
% plot(Mach_v, Ca_b)
% xlabel('Mach')
% ylabel('Ca_b')
% 
% subplot(2, 2, 4)
% plot(Mach_v, Ca_f)
% xlabel('Mach')
% ylabel('Ca_f')
% 
% 
% 
% figure
% subplot(2, 2, 1)
% plot(Mach_v, Cn)
% grid on
% xlabel('Mach')
% ylabel('Cn')
% legend('alpha: ',  int2str(alpha_v))
% 
% 
% subplot(2, 2, 2)
% plot(Mach_v, CL)
% grid on
% xlabel('Mach')
% ylabel('CL')
% legend('alpha: ',  int2str(alpha_v))
% 
% subplot(2, 2, 3)
% plot(Mach_v, Ca)
% grid on
% xlabel('Mach')
% ylabel('Ca')
% legend('alpha: ',  int2str(alpha_v))
% 
% subplot(2, 2, 4)
% plot(Mach_v, CD)
% grid on
% xlabel('Mach')
% ylabel('CD')
% legend('alpha: ',  int2str(alpha_v))
% 
% figure
% plot(Mach_v, CD, LineWidth=1.75)
% set(gca, 'FontSize', 30)
% xlabel('Mach', FontSize=35)
% ylabel('Cd', FontSize=35)
% legend('AoA = 0°', 'AoA = 5°', 'AoA = 10°', 'AoA = 15°', 'AoA = 20°', FontSize=30)
% axis square
% 
% figure
% plot(Mach_v, CL, LineWidth=1.75)
% set(gca, 'FontSize', 30)
% xlabel('Mach', FontSize=35)
% ylabel('Cl', FontSize=35)
% legend('AoA = 0°', 'AoA = 5°', 'AoA = 10°', 'AoA = 15°', 'AoA = 20°', FontSize=30)
% axis square

end

%% CDN
function Cdn = CDN(Mach, Mn_sample, cd_n_sample)
    % Interpolates the data extracted from the curve
    % Mn_sample: Mach values of the extracted points
    % cd_n_sample: Cdn values of the exctracted points
    
    C_p_stag = 1.8; % for Mach >= 4
    
    interp_fun = @(Mn) interp1(Mn_sample, cd_n_sample, Mn, 'pchip', 'extrap'); 
 
    
    if Mach < 4
        Cdn = interp_fun(Mach);
    else
        Cdn = (2/3) * C_p_stag;
    end
end
%% CN:
function [Cn,Cn_Cn0_sb,Cn_Cn0_Newt, Cdn, Cn_s] =  CN (a, b, x, ln, d_nose, A_b, A_r, A_p, M, alpha, h, phi, data, A_w)

% This function computes preliminary aerodynamic normal coefficient 
% according to the concept of component build-up
% inputs: 
% a, b = semimajor and semiminor axis  of cross section of the body [m]
% x = the coordinate at which cross section (a, b) are provided 
% x(1) = nosetip, x(end) = length of the body 
% d_nose = nose diameter 
% ln = length of the nose
% M, alpha, h = Mach, angle of attack, altitude.
% phi = cylindrical body orientation (0 or 90)
% 
% A_b --> area della base
% A_r --> area di riferimento (ref)
% 

eta = 0.7;      % varia in base al finesse ratio !!

if h > 84000
    Cn = 0;
end

alpha = deg2rad(alpha); % [rad]
l = x(end);
y = a;
R = y; % radius

fn_nose = ln/d_nose;

% deal with AoA
if alpha <= deg2rad(90) && alpha >= 0
    alpha = + alpha;
elseif alpha <= deg2rad(180) && alpha >= deg2rad(90)
    alpha = deg2rad(180) - alpha;
end

Cn_Cn0_sb = a/b*cos(deg2rad(phi))^2 + b/a*sin(deg2rad(phi))^2; % Slender-body experimental ratio

% if we want to consider a winged-body we can directly compute the
% following values by adding a term in the equation
if a == b
    Cn_Cn0_Newt = 1; 
else 
    if phi == 0
        Cn_Cn0_Newt = (3/2) * sqrt(a/b) * ((-b^2/a^2) / (1-(b^2/a^2)^(3/2)) * log((a/b) * (1+sqrt(1-(b^2/a^2))))+1/(1-(b^2/a^2)));
    elseif phi == 90
        Cn_Cn0_Newt = 3/2*sqrt(b/a)*(a^2/b^2/((a^2/b^2)-1)^(3/2)*atan(sqrt(a^2/b^2-1))-1/((a^2/b^2)-1)); 
    end
end

% % Excel file loading
% data = xlsread('Dataset Cdn.xlsx', 'Default Dataset');
% Estrazione delle coordinate x e y
Mn_sample = [0;data(:, 1);4]; % first column, Mach number 
Cd_n_sample = [1.2;data(:, 2);1.2863]; % Second column, Cdn

for i = 1:length(M)
    Cdn(i) = CDN(M(i), Mn_sample, Cd_n_sample); % Cross flow drag coefficient
end

% Body normal coefficient
Cn_b = A_b/A_r*sin(2*alpha)*cos(alpha/2)*Cn_Cn0_sb + eta*Cdn*A_p/A_r*sin(alpha)^2*Cn_Cn0_Newt; 

Cn = Cn_b;

A = x(end) / max(a); % body aspect ratio
% Tails
if nargin == 14
   if M^2 >= 1+(8/(pi*A))^2
       Cn_s = abs((4*abs(sin(alpha)*cos(alpha))/(M^2-1)^(1/2)+ 2*sin(alpha)^2)*(A_w/A_r));
   elseif M^2 < 1+(8/(pi*A))^2
       Cn_s = abs(((pi*A/2)*abs(sin(alpha)*cos(alpha))+2*sin(alpha)^2)*(A_w/A_r));
   end
   Cn = Cn_b + Cn_s;
end


end
%% drag_estimation:
function [Ca, Ca_w, Ca_f, Ca_b] = drag_estimation( a, b, x, ln, dn, nose_type, M, alpha, h, A_w)

% This function uses the method of component buildup
% to calculate the aerodynamic coefficients
% 
% a,b semimajor and semiminor axis of body
% 
% x scorre lungo l'asse
% x(1) --> nosetip
% x(n) --> total lenght
% 
% dnose
% lnose
% Mach
% alpha
% h
% fn --> finesse ratio
% 


if h > 84000
    Ca = 0;
end

alpha = deg2rad(alpha);     % deg to rad

l = x(end);
d = max(a);

y = a;

fn_nose = ln/dn;        % finesse ratio of nose cone



% Deal with AoA:

if alpha <= deg2rad(90) && alpha >= 0
    alpha = + alpha;
elseif alpha >= deg2rad(180) && alpha <= deg2rad(360)
    alpha = deg2rad(180);
end

if alpha <= deg2rad(90)
    theta = atan(0.5 / fn_nose);
elseif alpha > deg2rad(90)
    theta = pi/2;
end

switch nose_type
    case 'C'        % conical nose
        if M <= 1
            Ca_w = 0.8 * sin(theta)^2;
        elseif M > 1
            beta = sqrt(M^2-1);
            Ca_w = (4 * sin(theta)^2 * (2.5 + 8 * beta * sin(theta))) / (1 + 16 * beta * sin(theta));
        end

    case 'TO'       % tangent ogive nose
        if M <= 1
            Ca_w = 0.8 * sin(theta)^2;
        elseif M > 1
            Ca_w = wavedragogive(dn, ln, M);     % implementation of formula at pg 49
        end
end



% BASE PREASSURE:

gamma = 1.4;

% slide 57

if M >= 1

    Cp_b = 2/(gamma * M^2) * ( (2/(gamma+1))^1.4 * (1/M)^2.8 * (2 * gamma * M^2 - gamma + 1)/(gamma + 1) - 1 );

    Ca_b = -Cp_b;

elseif M < 1

    Ca_b = 0.12 + 0.13 * M^2;

end



% SKIN FRICTION:

% % Open rocket technical documentation ( consult )
% Cf_i = 1.48 * 10^-2;    % Re < 10^4
% % Cf_i = 0.0032 * (R_s/l)^2;        10^4<Re<10^5
% if M >= 1
%     C_sf = Cf_i/((1+0.144*M^2)^0.65);
% elseif M < 1
%     C_sf = Cf_i/(1+0.08*M^2);
% end
%%%%%%%%% % VAN DRIEST II: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if M == 0
    error('Mach = 0, non è possibile calcolare Cd_friction');
end
gamma = 1.4;
[T,a,P,rho,nu] = atmosisa(h);
V = M * a;                   % freestream velocity (m/s)

% Reynolds:
Re_x = V * d / nu;        % Reynolds number based on length

% Variables:
Tw_Taw = 1 + 0.9 * (gamma - 1)/2 * M.^2            % Wall-to-adiabatic temperature ratio
r = 0.88;                     % Recovery factor for supersonic flow
A = ( ( (gamma-1) * M.^2 ) ./ ( 2 .* Tw_Taw ) ).^0.5;
B = ( 1 + (gamma-1)/2 * M.^2 ) ./ ( Tw_Taw ) - 1;
C1 = ( 2*A.^2 - B ) ./ sqrt( B.^2 + 4*A.^2 );
C2 = B ./ sqrt( B.^2 + 4*A.^2 );

Cf = [];
for i = 1:length(M)
    f =@(x) log10(Re_x(i) .* x) - (1+2*r)/2 .* log10(Tw_Taw(i)) - ( 0.242 .* (asin(C1(i)) + asin(C2(i))) ) ./ ( A(i) .* x.^0.5 .* Tw_Taw(i).^0.5 );
    Cf = [Cf, fsolve(f, 0.01)];
end

% Mean turbulent skin friction global:
A_wet_body = (l-ln) * d * pi + (ln * dn * pi)/2 + A_w * 4;       % hp. cylindrical + 4 alette
A_ref = pi * d^2 / 4;
Cf_tot = Cf * A_wet_body/A_ref;


Ca_f = Cf_tot * (1 + 0.5 / (l/d) * A_wet_body) / A_ref;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% SUM:

Ca = Ca_w + Ca_f + Ca_b;
Ca = Ca * 1.1;      % there might be other types of drag sources (parasitic drag)

Ca = Ca*cos(alpha)^2;
end
%% wavedragogive:
function cdw = wavedragogive(diameter, lenght, Mach)

if Mach > 3.5
    Mach = 3.5;
end


sigma = 2 * (180/pi) *atan(diameter/2/lenght);
p =  (0.083  + 0.096/Mach^2) * (sigma/10)^1.69;
lod2 = (lenght/diameter)^2;
cdw = p * ( 1 - (196*lod2-16) / (14 * lod2 * (Mach+18)));
end
%% Body to Wind:
function [CL, CD] = body2wind(Cn, Ca, alpha_v)

CL = [];
CD = [];

CL = Cn .* cos(deg2rad(alpha_v)) - Ca .* sin(deg2rad(alpha_v));
CD = Cn .* sin(deg2rad(alpha_v)) + Ca .* cos(deg2rad(alpha_v));
end
%% Excel file loading:
function data = excel_load

% Excel file loading:
data = xlsread('Dataset Cdn.xlsx', 'Default Dataset');

end
