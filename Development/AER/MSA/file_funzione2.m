function [CD, CL_NKP, l_Cp_results, Cm] = file_funzione2(Mach_v, alpha_v, h, Xcg)

%
% La funzione calcola Cd e Cl prendendo in input: alpha, Mach, quota, e la geometria sotto, 
% superficie dell'area portante, posizione del centri di gravità (Xcg).
% E' i lementato Van Driest II per il calcolo della resistenza a frizione del corpo
% l'ala è modellata come 4 alette (90deg una dall'altra) immerse nel flusso
% INOLTRE:
% --> il calcolo usa una tecnica più corretta per la stima della portanza,
% --> calcolo della posizione del centro di pressione
% --> calcolo del momento aerodinamico
% 
% OSS: Xcg rispetto alla punta del razzo
% 

%% INPUTS:
geometry_funzione2
a_stall = 12+(6*cr_t/(2*s_t));
%% Calculation AXIAL COEFF.:

Ca = [];
Ca_w = [];
Ca_f = [];
Ca_b = [];


for i = 1:length(Mach_v)
    %Mach = Mach_v(i);

    for j = 1:length(alpha_v)

        %alpha = alpha_v(j);       

        [Ca(i,j), Ca_w(i,j), Ca_f(i,j), Ca_b(i,j)] = drag_estimation( a, b, x, ln, dn, nose_type, Mach_v(i), alpha_v(j), h, S_surface);

    end
end




%% Calculation NORMAL COEFF.:

Cn = [];

for i = 1:length(Mach_v)
    %Mach = Mach_v(i);

    for j = 1:length(alpha_v)

        %alpha = alpha_v(j);       

        [Cn(i,j),Cn_Cn0_sb(i,j),Cn_Cn0_Newt(i,j), Cdn(i,j),Cn_s(i,j)] =  CN (max(a), max(b), x, ln, dn, A_b, A_r, A_p, Mach_v(i), alpha_v(j), h, phi, data,S_surface);

    end
end

%% NKP

% Parametri per la convergenza
toll = 1e-2;
max_iter = 100;

for i = 1:length(Mach_v)
    for j = 1:length(alpha_v)
        % Inizializza variabili per il ciclo while
        l_Cp_old = 7;  % Ipotesi iniziale per il centro di pressione
        delta = inf;   % Delta iniziale
        iter = 0;

        while delta > toll && iter < max_iter
            iter = iter + 1;

            % Calcolo V_inf 
            a = 340; % velocità del suono a livello del mare [m/s]
            V_inf(i,j) = Mach_v(i)*a;


            % Calcolo del coefficiente di portanza
            x_hat_t = l_t - l_Cp_old;
             [CL_NKP(i,j),K_BT(i,j), K_TB(i,j),K_BW(i,j),K_WB(i,j),K_N(i,j),Cl_N(i,j),Cl_WB(i,j),Cl_BW(i,j),Cl_BT(i,j),Cl_TB(i,j),Cl_TV(i,j)] = Cl_NKP(alpha_v(j), Mach_v(i), h, r_N, S_w, s_w, r, cr_w, r_w, S_t,s_t,...
                           lambda_t, A_w, A_t,V_inf(i,j),r_t, l_w,l_t, x_hat_t,c_t, A_p);
           % Calcolo del lift totale
           rho = 1.225; % densità dell'aria a livello del mare
           L_tot(i,j) = 1/2*rho*V_inf(i,j)^2*S*CL_NKP(i,j);
           q_inf(i,j) = 1/2*rho*V_inf(i,j)^2;

            % Calcolo del centro di pressione
            [l_Cp_new] = l_CP(Mach_v(i), alpha_v(j), l_s,V_S,r_N,l_w,r, s_w,cr_w, q_inf(i,j), m_w, d, K_BW(i,j), K_WB(i,j), K_N(i,j), L_tot(i,j), s_t, l_t,...
                          cr_t,m_t, K_BT(i,j), K_TB(i,j), K_N(i,j), L_tot(i,j),Cl_N(i,j),Cl_WB(i,j), Cl_BW(i,j), Cl_BT(i,j), Cl_TB(i,j), Cl_TV(i,j));

            % Aggiorna delta
            delta = abs(l_Cp_new - l_Cp_old);

            % Aggiorna Cp_old
            l_Cp_old = l_Cp_new;
        end

        % % Salva i risultati finali
        % Cl_results(i, j) = Cl_NKP;
         l_Cp_results(i, j) = l_Cp_new;
    end
end


%% CL and CD conversion: (Jorgensen + Van Driest II)

[CL, CD] = body2wind(Cn, Ca, alpha_v);



%% Calcolo del coefficiente di momento aerodinamico:

d = 2 * r;      % diametro medio del body ( !! RIVEDERE !! )
Cm = (Xcg - l_Cp_results) .* Cn ./ d;



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

% figure
% plot(Mach_v, CL, LineWidth=1.75)
% set(gca, 'FontSize', 30)
% xlabel('Mach', FontSize=35)
% ylabel('Cl', FontSize=35)
% legend('AoA = 0°', 'AoA = 5°', 'AoA = 10°', 'AoA = 15°', 'AoA = 20°', FontSize=30)
% axis square
% hold on
% plot(Mach_v, CL_NKP)
% % xlabel('Mach')
% % ylabel('Cl_{NKP}')
% % legend('alpha: ',  int2str(alpha_v))
% 
% 
% figure
% plot(Mach_v, l_Cp_results, LineWidth=1.75)
% set(gca, 'FontSize', 30)
% xlabel('Mach', FontSize=35)
% ylabel('x_{Cp}', FontSize=35)
% legend('AoA = 0°', 'AoA = 5°', 'AoA = 10°', 'AoA = 15°', 'AoA = 20°', FontSize=30)
% % axis square


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
% A_w = surface area [m^2]
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
  Cn = (Cn_b*A_p + Cn_s*A_w)/(A_p+A_w);
end


end
%% drag_estimation:
function [Ca, Ca_w, Ca_f, Ca_b] = drag_estimation(a, b, x, ln, dn, nose_type, M, alpha, h, S_surface)

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
A_wet_body = (l-ln) * d * pi + (ln * dn * pi)/2 + S_surface * 2;       % hp. cylindrical + 2 * superficie bagnata
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

%% NKP function 
function [CL_NKP,K_BT, K_TB,K_BW,K_WB,K_N,Cl_N,Cl_WB,Cl_BW,Cl_BT,Cl_TB,Cl_TV,Cl_a_w, Cl_a_t] = Cl_NKP(alpha, M, h, r_N, S_w, s_w, r, cr_w, r_w, S_t,s_t,...
                           lambda_t, A_w, A_t,V_inf,r_t, l_w,l_t, x_hat_t,c_t,A_p)

% INPUT:
% alpha = angle of attack [deg]
% M = Mach number
% h = altitude [m]
% r_N = body radius at shoulder of nose [m]
% r = body radius [m]
% s = maximum semispsan of wing or tail in combination with body [m]
% S_w = wing-alone area [m^2]
% betaA = wing-alone or tail-alone effective aspect ratio
% lambda = taper ratio (0 = triangular wings, 1 = rectangular wings, 1/2 =
%                       trapezoidal wings)
% m = cotangent of leading edge sweep angle
% cr = chord at wing-body juncture or tail-body juncture [m]
% A_t = tail-alone aspect ratio
% A_w = wing-alone aspect ratio
% r_w = body radius at wing [m]
% V_inf = velocity [m/s]
% c_w = wing chord at midsection [m]
% Cl_w = lift coefficient based on the wing alone area
% r_t = body radius at tail [m]
% l_w = distance from most forward point of body to intersection of wing
%       leading edge and body [m]
% l_t = distance from most forward point of body to intersection of tail
%       leading edge and body [m]
% x_hat = distance to center of pressure measured from intersection of wing
%         leading edge and body for wing quantities and from intersection
%         of tail leading edge and body for tail quantities [m]
% c_t = tail chord at midsection [m]
% Cl_a_w = wing lift-curve slope for angle of attack, per radian
% Cl_a_t = tail lift-curve slope for angle of attack, per radian
% A_p = planform area [m^2]

if abs(alpha) < 1e-6  % Soglia di tolleranza
    alpha = 0.1;
end

d = 2*r;  % d = body diameter [m]

if h > 84000
    CL_NKP = 0;
end

alpha = deg2rad(alpha); % [rad]

% deal with AoA
if alpha <= deg2rad(90) && alpha >= 0
    alpha = + alpha;
elseif alpha <= deg2rad(180) && alpha >= deg2rad(90)
    alpha = deg2rad(180) - alpha;
end


% Cl_alpha

[Cl_a_w, Cl_a_t] = CL_A(M, A_w, A_t);

% Definition of Beta

if M == 1
        beta = sqrt(abs(0.9^2-1));
    else
        beta = sqrt(abs(M^2-1));
end


if s_w == 0 && s_t ~= 0  % BODY + TAIL (no wing)

    % LIFT ON BODY NOSE
    
    K_N = 2*pi*r_N^2/(S_t*Cl_a_t);
    Cl_N = K_N * Cl_a_t * alpha;

    % LIFT ON WING IN PRESENCE OF BODY
    K_WB = 0;
    Cl_WB = 0;

    % LIFT ON BODY DUE TO WING
    K_BW = 0;
    Cl_BW = 0;

    % LIFT ON TAIL IN PRESENCE OF BODY
    if r/s_t == 0 % the combination is all tail
        K_TB = 1;
    elseif r/s_t >= 1 % there is a very small exposed tail
        K_TB = 2;
    else
        K_TB = 2/pi*((1+r^4/s_t^4)*(1/2*atan(1/2*(s_t/r-r/s_t))+pi/4)-r^2/s_t^2*((s_t/r-r/s_t)+2*atan(r/s_t)))/(1-r/s_t)^2;
    end
    
    Cl_TB = K_TB * Cl_a_t * alpha;

    % LIFT ON BODY DUE TO TAIL
    
    if r/s_t == 0 % the combination is all tail
        K_BT = 0;
    elseif r/s_t >= 1 % there is a very small exposed tail
        K_BT = 2;
    else
        K_BT = ((1-r^2/s_t^2)^2 - 2/pi * ((1+r^4/s_t^4)*(1/2*atan(1/2*(s_t/r-r/s_t))+pi/4)-...
               r^2/s_t^2*((s_t/r-r/s_t)+2*atan(r/s_t))))/(1-r/s_t)^2;
    end
    
    Cl_BT = K_BT * Cl_a_t * alpha;

    
    % LIFT ON TAIL SECTIONS DUE TO WING VORTICES
    Cl_TV = 0;


    % TOTAL LIFT COEFFICIENT
   CL_NKP = (Cl_N*A_p+Cl_WB*S_w+Cl_BW*A_p+Cl_TB*S_t+Cl_BT*A_p+Cl_TV*S_t)/(S_w+S_t+A_p);

elseif s_t == 0 && s_w ~= 0  %  BODY + WING (no tail)
    
    % LIFT ON BODY NOSE
    
    K_N = 2*pi*r_N^2/(S_w*Cl_a_w);
    Cl_N = K_N * Cl_a_w * alpha;


    % LIFT ON WING IN PRESENCE OF BODY

    if r/s_w == 0 % the combination is all wing
        K_WB = 1;
    elseif r/s_w >= 1 % there is a very small exposed wing
        K_WB = 2;
    else
        K_WB = 2/pi*((1+r^4/s_w^4)*(1/2*atan(1/2*(s_w/r-r/s_w))+pi/4)-r^2/s_w^2*((s_w/r-r/s_w)+2*atan(r/s_w)))/(1-r/s_w)^2;
    end
    
    Cl_WB = K_WB * Cl_a_w * alpha;
    
    
    % LIFT ON BODY DUE TO WING
    
    if r/s_w == 0 % the combination is all wing
        K_BW = 0;
    elseif r/s_w >= 1 % there is a very small exposed wing
        K_BW = 2;
    else
        K_BW = ((1-r^2/s_w^2)^2 - 2/pi * ((1+r^4/s_w^4)*(1/2*atan(1/2*(s_w/r-r/s_w))+pi/4)-...
                r^2/s_w^2*((s_w/r-r/s_w)+2*atan(r/s_w))))/(1-r/s_w)^2;
    end

    Cl_BW = K_BW * Cl_a_w * alpha;

    % LIFT ON TAIL IN PRESENCE OF BODY (NEGLECTING WING VORTCES)
    K_TB = 0;
    Cl_TB = 0;

    % LIFT ON BODY DUE TAIL(NEGLECTING WING VORTICES)
    K_BT = 0;
    Cl_BT = 0;

    % LIFT ON TAIL SECTIONS DUE TO WING VORTICES
    Cl_TV = 0;


    % TOTAL LIFT COEFFICIENT
    CL_NKP = (Cl_N*A_p+Cl_WB*S_w+Cl_BW*A_p+Cl_TB*S_t+Cl_BT*A_p+Cl_TV*S_t)/(S_w+S_t+A_p);
elseif s_w ~= 0 && s_t ~= 0  % BODY + WING + TAIL

    % LIFT ON BODY NOSE
    
    K_N = 2*pi*r_N^2/(S_w*Cl_a_w);
    Cl_N = K_N * Cl_a_w * alpha;


    % LIFT ON WING IN PRESENCE OF BODY

    if r/s_w == 0 % the combination is all wing
        K_WB = 1;
    elseif r/s_w >= 1 % there is a very small exposed wing
        K_WB = 2;
    else
        K_WB = 2/pi*((1+r^4/s_w^4)*(1/2*atan(1/2*(s_w/r-r/s_w))+pi/4)-r^2/s_w^2*((s_w/r-r/s_w)+2*atan(r/s_w)))/(1-r/s_w)^2;
    end
    
    Cl_WB = K_WB * Cl_a_w * alpha;
    
    
    % LIFT ON BODY DUE TO WING
    
    if r/s_w == 0 % the combination is all wing
        K_BW = 0;
    elseif r/s_w >= 1 % there is a very small exposed wing
        K_BW = 2;
    else
        K_BW = ((1-r^2/s_w^2)^2 - 2/pi * ((1+r^4/s_w^4)*(1/2*atan(1/2*(s_w/r-r/s_w))+pi/4)-...
               r^2/s_w^2*((s_w/r-r/s_w)+2*atan(r/s_w))))/(1-r/s_w)^2;
    end
    
    Cl_BW = K_BW * Cl_a_w * alpha;
    
    
    % LIFT ON TAIL IN PRESENCE OF BODY (NEGLECTING WING VORTCES)
    
    if r/s_t == 0 % the combination is all tail
        K_TB = 1;
    elseif r/s_t >= 1 % there is a very small exposed tail
        K_TB = 2;
    else
        K_TB = 2/pi*((1+r^4/s_t^4)*(1/2*atan(1/2*(s_t/r-r/s_t))+pi/4)-r^2/s_t^2*((s_t/r-r/s_t)+2*atan(r/s_t)))/(1-r/s_t)^2;
    end
    
    Cl_TB = K_TB * Cl_a_t * alpha * S_t/S_w;
    
    
    
    % LIFT ON BODY DUE TAIL(NEGLECTING WING VORTICES)
    
        if r/s_t == 0 % the combination is all tail
            K_BT = 0;
        elseif r/s_t >= 1 % there is a very small exposed tail
            K_BT = 2;
        else
           K_BT = ((1-r^2/s_t^2)^2 - 2/pi * ((1+r^4/s_t^4)*(1/2*atan(1/2*(s_t/r-r/s_t))+pi/4)-...
               r^2/s_t^2*((s_t/r-r/s_t)+2*atan(r/s_t))))/(1-r/s_t)^2;
        end
    
    Cl_BT = K_BT * Cl_a_t * alpha * S_t/S_w;
    
    
    % LIFT ON TAIL SECTIONS DUE TO WING VORTICES
    
    cl = 0.1; % giustificare con i paper
    Cl_w = Cl_a_w*alpha;
    Cl_t = Cl_a_t*alpha;
    f_r_s_r_w = (pi/4-pi/4*(r/s_w)^2-(r/s_w)+(1+(r/s_w)^2)^2/(2*(1-(r/s_w)^2))*asin((1-(r/s_w)^2)/(1+(r/s_w)^2)))/...
                (2*(1-(r/s_w)));
    f_w =f_r_s_r_w*(s_w-r)+r; 
    f_t = Cl_t*S_t/(2*cl*c_t);
    Gamma_m = V_inf*K_WB*alpha*Cl_a_w*S_w/(4*(f_w-r_w));
    h_t = (l_t+x_hat_t-l_w-cr_w)*sin(alpha);
    f_i = f_w*r_t^2/(f_w^2+h_t^2);
    f_i_ = -f_i;
    f_w_ = -f_w;
    h_i = h_t*r_t^2/(f_w^2+h_t^2);
    [L1] = L(lambda_t, r_t, s_t, f_w, h_t);
    [L2] = L(lambda_t, r_t, s_t, f_w_, h_t);
    [L3] = L(lambda_t, r_t, s_t, f_i, h_i);
    [L4] = L(lambda_t, r_t, s_t, f_i_, h_i);
    
    i = 2/(1+lambda_t)*(L1-L2-L3+L4);
    
    
    Cl_TV = Cl_a_w*Cl_a_t*K_WB*alpha*i*(s_t-r_t)/(2*pi*A_t*(f_w-r_w));
    
    
    % TOTAL LIFT COEFFICIENT
    
    CL_NKP = (Cl_N*A_p+Cl_WB*S_w+Cl_BW*A_p+Cl_TB*S_t+Cl_BT*A_p+Cl_TV*S_t)/(S_w+S_t+A_p);
end


end


%% L function 

function [L] = L(lambda, r, s, f, h)

% lambda = taper ratio
% r = radius of body at tail[m]
% s = maximum semispsan of tail in combination with body [m]
% f = wing vortex semispan at wing trailing edge [m]
% h = height of wing vortex above body axis at tail centre of pressure [m]


L = (((s-lambda*r)-f*(1-lambda))/(2*(s-r))*log((h^2+(f-s)^2)/(h^2+(f-r)^2))-...
    (1-lambda)/(s-r)*((s-r)+h*atan((f-s)/h)-h*atan((f-r)/h)));

end


%% l_CP function 

function [l_CP] = l_CP(M, alpha, l_s,V_S,r_N,l_w,r, s_w,cr_w, q_inf, m_w, d, K_BW, K_WB, K_N_W, LW, s_t, l_t,...
    cr_t,m_t, K_BT, K_TB, K_N_T, LT,Cl_N,Cl_WB, Cl_BW, Cl_BT, Cl_TB, Cl_TV)

% INPUTS:
% M = mach number
% alpha = angle of attack
% l_s = distance from most forward point of body to shoulder of body  nose [m]
% V_S = volume of body nose up to shoulder [m^3]
% r_N = body radius at shoulder of nose [m]
% l_w = distance from most forward point of body to intersection of wing
%       leading edge and body [m]
% r = body radius [m]
% s = maximum semispan of wing or tail in combination with body [m]
% cr = chord at wing-body juncture or tail-body juncture [m]
% q_inf = free stream dynamic pressure [kg/m^2]
% m = cotangent of leading edge sweep angle
% d = body diameter [m]
% K_BW = ratio of lift component of body in presence of wing
% K_WB = ratio of lift component of wing in presence of body
% K_N_W = ratio of lift of body nose
% LW = lift force of the wing-body 
% l_t = distance from most forward point of body to intersection of tail
%       leading edge and body [m]
% K_BT = ratio of lift component of body in presence of tail
% K_TB = ratio of lift component of tail in presence of body
% K_N_T = ratio of lift of body nose
% LT = lift force of the tail-body 
% Cl_N = lift coefficient based on body nose
% Cl_WB = lift coefficient based on wing in presence of body
% Cl_BW = lift coefficient based on body in presence of wing
% Cl_BT = lift coefficient based on body in presence of tail
% Cl_TB = lift coefficient based on tail in presence of body
% Cl_TV = lift coefficient based on tail in presence of wing vortex

if abs(alpha) < 1e-6  % Soglia di tolleranza
    alpha = 0.1;
end


alpha = deg2rad(alpha);


% alpha = deg2rad(alpha); % [rad]

% deal with AoA
if alpha <= deg2rad(90) && alpha >= 0
    alpha = + alpha;
elseif alpha <= deg2rad(180) && alpha >= deg2rad(90)
    alpha = deg2rad(180) - alpha;
end


% Definition of Beta

if M == 1
    beta = sqrt(abs(0.9^2-1));
else
    beta = sqrt(abs(M^2-1));
end


% CENTRE OF PRESSURE OF BODY NOSE

l_N = l_s*(1-V_S/(pi*r_N^2*l_s));



if s_w == 0 && s_t ~= 0 % BODY + TAIL (no wing)
    
    % CENTRE OF PRESSURE OF WING IN PRESENCE OF BODY
    l_WB = 0;

    % CENTRE OF PRESSURE ON BODY DUE TO WING
    l_BW = 0;

    % CENTRE OF PRESSURE OF TAIL IN PRESENCE OF BODY

    x_cr_TB_a = 1/(1-r/s_t)*(2*(1/3+r^4/s_t^4)*atan(s_t/r)+2/3*r^3/s_t^3*log(((s_t^2+r^2)/(2*s_t^2))^2*s_t/r)-1/3*r^3/s_t^3*...
                (2*pi-1+s_t^2/r^2))/((1+r^2/s_t^2)^2*atan(s_t/r)-r^2/s_t^2*(pi+(s_t/r-r/s_t)))-r/s_t/(1-r/s_t);
    
    l_TB = l_t + cr_t*x_cr_TB_a*alpha;

    % CENTRE OF PRESSURE ON BODY DUE TO TAIL

    if M > 1
    
        M_BT = 4*q_inf*alpha*m_t/(3*pi*beta)*cr_t^3*(sqrt(1+2*beta*d/cr_t)*((2*m_t*beta+5)/(3*(m_t*beta+1)^2)+beta*d/cr_t/(3*...
            (m_t*beta+1))-(beta*d/cr_t)^2/(beta*m_t))+1/sqrt(m_t^2*beta^2-1)*((1+beta*d/cr_t)^3-(beta*d/cr_t)^3/(m_t^2*beta^2)-1/...
            (1+m_t*beta)^2)*acos((1+beta*d/cr_t*(m_t*beta+1))/(m_t*beta+beta*d/cr_t*(m_t*beta+1)))+(beta*d/cr_t)^3*1/(m_t^2*beta^2)*...
            acosh(1+cr_t/(beta*d))-((2*m_t*beta+5)/(3*(m_t*beta+1)^2))-(1-(1/(m_t*beta+1))^2)/sqrt(m_t^2*beta^2-1)*acos(1/(m_t*beta)));
    
    elseif M <= 1
    
        M_BT = 4*q_inf*alpha/(pi*beta^2)*cr_t^3*(sqrt(m_t^2*beta^2+m_t*beta*(m_t*beta+1)*beta*d/cr_t)/(9*m_t*beta*(m_t*beta+1)^3)*...
            ((8*m_t*beta+24)*m_t^2*beta^2+(14*m_t*beta+6)*(m_t*beta+1)*m_t*beta*beta*d/cr_t+3*(m_t*beta-3)*(m_t*beta+1)^2*(beta*d/cr_t)^2)-...
            (8*m_t*beta+24)*m_t^3*beta^3/(9*m_t*beta*(m_t*beta+1)^3)-(m_t*beta-3)/(3*m_t*beta)*(beta*d/cr_t)^3*acosh(sqrt((m_t*beta+(m_t*beta+1)*...
            beta*d/cr_t)/((m_t*beta+1)*beta*d/cr_t))));
    end
    
    K_CT = K_BT+K_TB+K_N_T;
    L_T = LT/K_CT; % lift of the tail alone
    
    x_cr_BT = M_BT/(K_BT*L_T*cr_t);
    
    l_BT = l_t+cr_t*x_cr_BT;

    % CENTRE OF PRESSURE OF TAIL DUE TO WING VORTICES

    l_TV = 0;


    
    % CENTRE OF PRESSURE FOR ENTIRE COMBINATION
    
    l_CP = (l_N*Cl_N+l_WB*Cl_WB+l_BW*Cl_BW+l_BT*Cl_BT+l_TB*Cl_TB+l_TV*Cl_TV)/(Cl_N+Cl_WB+Cl_BW+Cl_BT+Cl_TB+Cl_TV);


elseif s_t == 0 && s_w ~= 0 % BODY + WING (no tail)
    
    % CENTRE OF PRESSURE OF WING IN PRESENCE OF BODY

    x_cr_WB_a = 1/(1-r/s_w)*(2*(1/3+r^4/s_w^4)*atan(s_w/r)+2/3*r^3/s_w^3*log(((s_w^2+r^2)/(2*s_w^2))^2*s_w/r)-1/3*r^3/s_w^3*...
        (2*pi-1+s_w^2/r^2))/((1+r^2/s_w^2)^2*atan(s_w/r)-r^2/s_w^2*(pi+(s_w/r-r/s_w)))-r/s_w/(1-r/s_w);
    
    l_WB = l_w + cr_w*x_cr_WB_a*alpha;
    
    
    
    % CENTRE OF PRESSURE ON BODY DUE TO WING
      
    if M > 1
    
        M_BW = 4*q_inf*alpha*m_w/(3*pi*beta)*cr_w^3*(sqrt(1+2*beta*d/cr_w)*((2*m_w*beta+5)/(3*(m_w*beta+1)^2)+beta*d/cr_w/(3*...
            (m_w*beta+1))-(beta*d/cr_w)^2/(beta*m_w))+1/sqrt(m_w^2*beta^2-1)*((1+beta*d/cr_w)^3-(beta*d/cr_w)^3/(m_w^2*beta^2)-1/...
            (1+m_w*beta)^2)*acos((1+beta*d/cr_w*(m_w*beta+1))/(m_w*beta+beta*d/cr_w*(m_w*beta+1)))+(beta*d/cr_w)^3*1/(m_w^2*beta^2)*...
            acosh(1+cr_w/(beta*d))-((2*m_w*beta+5)/(3*(m_w*beta+1)^2))-(1-(1/(m_w*beta+1))^2)/sqrt(m_w^2*beta^2-1)*acos(1/(m_w*beta)));
    
    elseif M <= 1
    
        M_BW = 4*q_inf*alpha/(pi*beta^2)*cr_w^3*(sqrt(m_w^2*beta^2+m_w*beta*(m_w*beta+1)*beta*d/cr_w)/(9*m_w*beta*(m_w*beta+1)^3)*...
            ((8*m_w*beta+24)*m_w^2*beta^2+(14*m_w*beta+6)*(m_w*beta+1)*m_w*beta*beta*d/cr_w+3*(m_w*beta-3)*(m_w*beta+1)^2*(beta*d/cr_w)^2)-...
            (8*m_w*beta+24)*m_w^3*beta^3/(9*m_w*beta*(m_w*beta+1)^3)-(m_w*beta-3)/(3*m_w*beta)*(beta*d/cr_w)^3*acosh(sqrt((m_w*beta+(m_w*beta+1)*...
            beta*d/cr_w)/((m_w*beta+1)*beta*d/cr_w))));
    end
    
    K_CW = K_BW+K_WB+K_N_W;
    L_W = LW/K_CW; % lift of the wing alone
    
    x_cr_BW = M_BW/(K_BW*L_W*cr_w);
    
    l_BW = l_w+cr_w*x_cr_BW;

    % CENTRE OF PRESSURE OF TAIL IN PRESENCE OF BODY
    l_TB = 0;
    
    % CENTRE OF PRESSURE ON BODY DUE TO TAIL
    l_BT = 0;

    % CENTRE OF PRESSURE OF TAIL DUE TO WING VORTICES

    l_TV = l_TB;
    
    
    
    % CENTRE OF PRESSURE FOR ENTIRE COMBINATION
    
    l_CP = (l_N*Cl_N+l_WB*Cl_WB+l_BW*Cl_BW+l_BT*Cl_BT+l_TB*Cl_TB+l_TV*Cl_TV)/(Cl_N+Cl_WB+Cl_BW+Cl_BT+Cl_TB+Cl_TV);

elseif s_w ~= 0 && s_t ~= 0 % BODY + TAIL + WING

    % CENTRE OF PRESSURE OF WING IN PRESENCE OF BODY
    
    x_cr_WB_a = 1/(1-r/s_w)*(2*(1/3+r^4/s_w^4)*atan(s_w/r)+2/3*r^3/s_w^3*log(((s_w^2+r^2)/(2*s_w^2))^2*s_w/r)-1/3*r^3/s_w^3*...
        (2*pi-1+s_w^2/r^2))/((1+r^2/s_w^2)^2*atan(s_w/r)-r^2/s_w^2*(pi+(s_w/r-r/s_w)))-r/s_w/(1-r/s_w);
    
    l_WB = l_w + cr_w*x_cr_WB_a*alpha;
    
        
    % CENTRE OF PRESSURE ON BODY DUE TO WING
    
    if M > 1
    
        M_BW = 4*q_inf*alpha*m_w/(3*pi*beta)*cr_w^3*(sqrt(1+2*beta*d/cr_w)*((2*m_w*beta+5)/(3*(m_w*beta+1)^2)+beta*d/cr_w/(3*...
            (m_w*beta+1))-(beta*d/cr_w)^2/(beta*m_w))+1/sqrt(m_w^2*beta^2-1)*((1+beta*d/cr_w)^3-(beta*d/cr_w)^3/(m_w^2*beta^2)-1/...
            (1+m_w*beta)^2)*acos((1+beta*d/cr_w*(m_w*beta+1))/(m_w*beta+beta*d/cr_w*(m_w*beta+1)))+(beta*d/cr_w)^3*1/(m_w^2*beta^2)*...
            acosh(1+cr_w/(beta*d))-((2*m_w*beta+5)/(3*(m_w*beta+1)^2))-(1-(1/(m_w*beta+1))^2)/sqrt(m_w^2*beta^2-1)*acos(1/(m_w*beta)));
    
    elseif M <= 1
    
        M_BW = 4*q_inf*alpha/(pi*beta^2)*cr_w^3*(sqrt(m_w^2*beta^2+m_w*beta*(m_w*beta+1)*beta*d/cr_w)/(9*m_w*beta*(m_w*beta+1)^3)*...
            ((8*m_w*beta+24)*m_w^2*beta^2+(14*m_w*beta+6)*(m_w*beta+1)*m_w*beta*beta*d/cr_w+3*(m_w*beta-3)*(m_w*beta+1)^2*(beta*d/cr_w)^2)-...
            (8*m_w*beta+24)*m_w^3*beta^3/(9*m_w*beta*(m_w*beta+1)^3)-(m_w*beta-3)/(3*m_w*beta)*(beta*d/cr_w)^3*acosh(sqrt((m_w*beta+(m_w*beta+1)*...
            beta*d/cr_w)/((m_w*beta+1)*beta*d/cr_w))));
    end
    
    K_CW = K_BW+K_WB+K_N_W;
    L_W = LW/K_CW; % lift of the wing alone
    
    x_cr_BW = M_BW/(K_BW*L_W*cr_w);
    
    l_BW = l_w+cr_w*x_cr_BW;
    
    
    % CENTRE OF PRESSURE OF TAIL IN PRESENCE OF BODY
    
    x_cr_TB_a = 1/(1-r/s_t)*(2*(1/3+r^4/s_t^4)*atan(s_t/r)+2/3*r^3/s_t^3*log(((s_t^2+r^2)/(2*s_t^2))^2*s_t/r)-1/3*r^3/s_t^3*...
        (2*pi-1+s_t^2/r^2))/((1+r^2/s_t^2)^2*atan(s_t/r)-r^2/s_t^2*(pi+(s_t/r-r/s_t)))-r/s_t/(1-r/s_t);
    
    l_TB = l_t + cr_t*x_cr_TB_a*alpha;
    
    
    % CENTRE OF PRESSURE ON BODY DUE TO TAIL
    
    if M > 1
    
        M_BT = 4*q_inf*alpha*m_t/(3*pi*beta)*cr_t^3*(sqrt(1+2*beta*d/cr_t)*((2*m_t*beta+5)/(3*(m_t*beta+1)^2)+beta*d/cr_t/(3*...
            (m_t*beta+1))-(beta*d/cr_t)^2/(beta*m_t))+1/sqrt(m_t^2*beta^2-1)*((1+beta*d/cr_t)^3-(beta*d/cr_t)^3/(m_t^2*beta^2)-1/...
            (1+m_t*beta)^2)*acos((1+beta*d/cr_t*(m_t*beta+1))/(m_t*beta+beta*d/cr_t*(m_t*beta+1)))+(beta*d/cr_t)^3*1/(m_t^2*beta^2)*...
            acosh(1+cr_t/(beta*d))-((2*m_t*beta+5)/(3*(m_t*beta+1)^2))-(1-(1/(m_t*beta+1))^2)/sqrt(m_t^2*beta^2-1)*acos(1/(m_t*beta)));
    
    elseif M <= 1
    
        M_BT = 4*q_inf*alpha/(pi*beta^2)*cr_t^3*(sqrt(m_t^2*beta^2+m_t*beta*(m_t*beta+1)*beta*d/cr_t)/(9*m_t*beta*(m_t*beta+1)^3)*...
            ((8*m_t*beta+24)*m_t^2*beta^2+(14*m_t*beta+6)*(m_t*beta+1)*m_t*beta*beta*d/cr_t+3*(m_t*beta-3)*(m_t*beta+1)^2*(beta*d/cr_t)^2)-...
            (8*m_t*beta+24)*m_t^3*beta^3/(9*m_t*beta*(m_t*beta+1)^3)-(m_t*beta-3)/(3*m_t*beta)*(beta*d/cr_t)^3*acosh(sqrt((m_t*beta+(m_t*beta+1)*...
            beta*d/cr_t)/((m_t*beta+1)*beta*d/cr_t))));
    end
    
    K_CT = K_BT+K_TB+K_N_T;
    L_T = LT/K_CT; % lift of the tail alone
    
    x_cr_BT = M_BT/(K_BT*L_T*cr_t);
    
    l_BT = l_t+cr_t*x_cr_BT;
    
    % CENTRE OF PRESSURE OF TAIL DUE TO WING VORTICES
    
    l_TV = l_TB;
    
    % CENTRE OF PRESSURE FOR ENTIRE COMBINATION
    
    l_CP = (l_N*Cl_N+l_WB*Cl_WB+l_BW*Cl_BW+l_BT*Cl_BT+l_TB*Cl_TB+l_TV*Cl_TV)/(Cl_N+Cl_WB+Cl_BW+Cl_BT+Cl_TB+Cl_TV);

end


end


%% CL_A function 
function [Cl_a_w, Cl_a_t] = CL_A(M, A_w, A_t)

if A_w == 0 && A_t ~= 0
    Cl_a_w = 0;

     if M < 0.3
   
    Cl_a_t = 2*pi*0.4/(1+2*0.7/A_t);
    elseif M >= 0.3 && M < 1
       
        Cl_a_t_inc = 2*pi*0.4/(1+2*0.7/A_t);
       
        Cl_a_t = Cl_a_t_inc/sqrt(1-M^2); % Correzione di Prandtl-Glauert per effetti di compressibilità
    elseif M == 1
         
         Cl_a_t_inc = 2*pi*0.67/(1+2*0.8/A_t);
        
         Cl_a_t = Cl_a_t_inc/sqrt(1-0.92^2); 
    elseif M > 1 && M <= 1.1
        Cl_a_t = 4*0.85/sqrt(M^2-1);
       
    elseif M > 1.1 && M <= 1.2
        Cl_a_t = 4*0.95/sqrt(M^2-1);
        
    elseif M > 1.2 && M <= 1.3
        Cl_a_t = 4*1.00/sqrt(M^2-1);
        
    elseif M > 1.3 && M <= 1.4
        Cl_a_t = 4*1.05/sqrt(M^2-1);
        
    elseif M > 1.4 && M <= 1.5
        Cl_a_w = 4*1.1/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 1.5 && M <= 1.6
        Cl_a_w = 4*1.16/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 1.6 && M <= 1.7
        Cl_a_w = 4*1.16/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 1.7 && M <= 1.8
        Cl_a_w = 4*1.21/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 1.8 && M <= 1.9
        Cl_a_w = 4*1.25/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 1.9 && M <= 2
        Cl_a_w = 4*1.29/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2 && M <= 2.1
        Cl_a_w = 4*1.34/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2.1 && M <= 2.2
        Cl_a_w = 4*1.37/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2.2 && M <= 2.3
        Cl_a_w = 4*1.39/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2.3 && M <= 2.4
        Cl_a_w = 4*1.40/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2.4 && M <= 2.5
        Cl_a_w = 4*1.41/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2.5 && M <= 2.6
        Cl_a_w = 4*1.42/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2.6 && M <= 2.7
        Cl_a_w = 4*1.43/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2.7 && M <= 2.8
        Cl_a_w = 4*1.44/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2.8 && M <= 2.9
        Cl_a_w = 4*1.45/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2.9 && M <= 3
        Cl_a_w = 4*1.46/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3 && M <= 3.1
        Cl_a_w = 4*1.47/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3.1 && M <= 3.2
        Cl_a_w = 4*1.48/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3.2 && M <= 3.3
        Cl_a_w = 4*1.50/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3.3 && M <= 3.4
        Cl_a_w = 4*1.52/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3.4 && M <= 3.5
        Cl_a_w = 4*1.54/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3.5 && M <= 3.6
        Cl_a_w = 4*1.56/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3.6 && M <= 3.7
        Cl_a_w = 4*1.58/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3.7 && M <= 3.8
        Cl_a_w = 4*1.60/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3.8 && M <= 3.9
        Cl_a_w = 4*1.62/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3.9 && M <= 4
        Cl_a_w = 4*1.64/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4 && M <= 4.1
        Cl_a_w = 4*1.66/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4.1 && M <= 4.2
        Cl_a_w = 4*1.68/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4.2 && M <= 4.3
        Cl_a_w = 4*1.70/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4.3 && M <= 4.4
        Cl_a_w = 4*1.73/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4.4 && M <= 4.5
        Cl_a_w = 4*1.76/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4.5 && M <= 4.6
        Cl_a_w = 4*1.79/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4.6 && M <= 4.7
        Cl_a_w = 4*1.82/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4.7 && M <= 4.8
        Cl_a_w = 4*1.85/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4.8 && M <= 4.9
        Cl_a_w = 4*1.88/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4.9 && M <= 5
        Cl_a_w = 4*1.91/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5 && M <= 5.1
        Cl_a_w = 4*1.94/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5.1 && M <= 5.2
        Cl_a_w = 4*1.97/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5.2 && M <= 5.3
        Cl_a_w = 4*2.00/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5.3 && M <= 5.4
        Cl_a_w = 4*2.03/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5.4 && M <= 5.5
        Cl_a_w = 4*2.06/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5.5 && M <= 5.6
        Cl_a_w = 4*2.09/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5.6 && M <= 5.7
        Cl_a_w = 4*2.12/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5.7 && M <= 5.8
        Cl_a_w = 4*2.15/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5.8 && M <= 5.9
        Cl_a_w = 4*2.18/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5.9 && M <= 6
        Cl_a_w = 4*2.21/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6 && M <= 6.1
        Cl_a_w = 4*2.24/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6.1 && M <= 6.2
        Cl_a_w = 4*2.27/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6.2 && M <= 6.3
        Cl_a_w = 4*2.30/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6.3 && M <= 6.4
        Cl_a_w = 4*2.33/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6.4 && M <= 6.5
        Cl_a_w = 4*2.36/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6.5 && M <= 6.6
        Cl_a_w = 4*2.39/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6.6 && M <= 6.7
        Cl_a_w = 4*2.42/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6.7 && M <= 6.8
        Cl_a_w = 4*2.45/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6.8 && M <= 6.9
        Cl_a_w = 4*2.48/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6.9 && M <= 7
        Cl_a_w = 4*2.51/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 7 && M <= 7.1
        Cl_a_w = 4*2.54/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 7.1 && M <= 7.2
        Cl_a_w = 4*2.57/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 7.2 && M <= 7.3
        Cl_a_w = 4*2.60/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 7.3 && M <= 7.4
        Cl_a_w = 4*2.63/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 7.4 
        Cl_a_w = 4*2.66/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    end

elseif A_w ~= 0 && A_t == 0
    Cl_a_t = 0;

    if M < 0.3
    Cl_a_t = 2*pi*0.4/(1+2*0.7/A_t);

    elseif M >= 0.3 && M < 1
        Cl_a_w_inc = 2*pi*0.4/(1+2*0.7/A_t);
        Cl_a_w = Cl_a_w_inc/sqrt(1-M^2); % Correzione di Prandtl-Glauert per effetti di compressibilità
    
    elseif M == 1
        Cl_a_w_inc = 2*pi*0.4/(1+2*0.8/A_t);
        Cl_a_w = Cl_a_w_inc/sqrt(1-0.92^2); 
    
    elseif M > 1 && M <= 3
        Cl_a_w = 4/sqrt(M^2-1);
        
    elseif M > 3 && M <= 3.1
        Cl_a_w = 4*1.02/sqrt(M^2-1);
        
    elseif M > 3.1 && M <= 3.2
        Cl_a_w = 4*1.04/sqrt(M^2-1);
        
    elseif M > 3.2 && M <= 3.3
        Cl_a_w = 4*1.06/sqrt(M^2-1);
        
    elseif M > 3.3 && M <= 3.4
        Cl_a_w = 4*1.08/sqrt(M^2-1);
        
    elseif M > 3.4 && M <= 3.5
        Cl_a_w = 4*1.10/sqrt(M^2-1);
        
    elseif M > 3.5 && M <= 3.6
        Cl_a_w = 4*1.12/sqrt(M^2-1);
        
    elseif M > 3.6 && M <= 3.7
        Cl_a_w = 4*1.14/sqrt(M^2-1);
        
    elseif M > 3.7 && M <= 3.8
        Cl_a_w = 4*1.16/sqrt(M^2-1);
        
    elseif M > 3.8 && M <= 3.9
        Cl_a_w = 4*1.18/sqrt(M^2-1);
        
    elseif M > 3.9 && M <= 4
        Cl_a_w = 4*1.20/sqrt(M^2-1);
        
    elseif M > 4 && M <= 4.1
        Cl_a_w = 4*1.22/sqrt(M^2-1);
        
    elseif M > 4.1 && M <= 4.2
        Cl_a_w = 4*1.24/sqrt(M^2-1);
        
    elseif M > 4.2 && M <= 4.3
        Cl_a_w = 4*1.26/sqrt(M^2-1);
        
    elseif M > 4.3 && M <= 4.4
        Cl_a_w = 4*1.28/sqrt(M^2-1);
        
    elseif M > 4.4 && M <= 4.5
        Cl_a_w = 4*1.30/sqrt(M^2-1);
        
    elseif M > 4.5 && M <= 4.6
        Cl_a_w = 4*1.32/sqrt(M^2-1);
        
    elseif M > 4.6 && M <= 4.7
        Cl_a_w = 4*1.34/sqrt(M^2-1);
        
    elseif M > 4.7 && M <= 4.8
        Cl_a_w = 4*1.36/sqrt(M^2-1);
        
    elseif M > 4.8 && M <= 4.9
        Cl_a_w = 4*1.38/sqrt(M^2-1);
        
    elseif M > 4.9 && M <= 5
        Cl_a_w = 4*1.40/sqrt(M^2-1);
        
    elseif M > 5 && M <= 5.1
        Cl_a_w = 4*1.42/sqrt(M^2-1);
        
    elseif M > 5.1 && M <= 5.2
        Cl_a_w = 4*1.44/sqrt(M^2-1);
        
    elseif M > 5.2 && M <= 5.3
        Cl_a_w = 4*1.46/sqrt(M^2-1);
        
    elseif M > 5.3 && M <= 5.4
        Cl_a_w = 4*1.48/sqrt(M^2-1);
        
    elseif M > 5.4 && M <= 5.5
        Cl_a_w = 4*1.50/sqrt(M^2-1);
        
    elseif M > 5.5 && M <= 5.6
        Cl_a_w = 4*1.52/sqrt(M^2-1);
        
    elseif M > 5.6 && M <= 5.7
        Cl_a_w = 4*1.54/sqrt(M^2-1);
        
    elseif M > 5.7 && M <= 5.8
        Cl_a_w = 4*1.56/sqrt(M^2-1);
        
    elseif M > 5.8 && M <= 5.9
        Cl_a_w = 4*1.58/sqrt(M^2-1);
        
    elseif M > 5.9 && M <= 6
        Cl_a_w = 4*1.60/sqrt(M^2-1);
        
    elseif M > 6 && M <= 6.1
        Cl_a_w = 4*1.62/sqrt(M^2-1);
        
    elseif M > 6.1 && M <= 6.2
        Cl_a_w = 4*1.64/sqrt(M^2-1);
        
    elseif M > 6.2 && M <= 6.3
        Cl_a_w = 4*1.66/sqrt(M^2-1);
        
    elseif M > 6.3 && M <= 6.4
        Cl_a_w = 4*1.68/sqrt(M^2-1);
        
    elseif M > 6.4 && M <= 6.5
        Cl_a_w = 4*1.70/sqrt(M^2-1);
        
    elseif M > 6.5 && M <= 6.6
        Cl_a_w = 4*1.72/sqrt(M^2-1);
       
    elseif M > 6.6 && M <= 6.7
        Cl_a_w = 4*1.74/sqrt(M^2-1);
        
    elseif M > 6.7 && M <= 6.8
        Cl_a_w = 4*1.76/sqrt(M^2-1);
        
    elseif M > 6.8 && M <= 6.9
        Cl_a_w = 4*1.78/sqrt(M^2-1);
        
    elseif M > 6.9 && M <= 7
        Cl_a_w = 4*1.80/sqrt(M^2-1);
        
    elseif M > 7 && M <= 7.1
        Cl_a_w = 4*1.82/sqrt(M^2-1);
        
    elseif M > 7.1 && M <= 7.2
        Cl_a_w = 4*1.84/sqrt(M^2-1);
        
    elseif M > 7.2 && M <= 7.3
        Cl_a_w = 4*1.86/sqrt(M^2-1);
        
    elseif M > 7.3 && M <= 7.4
        Cl_a_w = 4*1.88/sqrt(M^2-1);
        
    elseif M > 7.4 
        Cl_a_w = 4*1.90/sqrt(M^2-1);
        
    end

elseif A_w ~= 0 && A_t ~= 0

    if M < 0.3
    Cl_a_w = 2*pi*0.4/(1+2*0.7/A_w);
    Cl_a_t = 2*pi*0.4/(1+2*0.7/A_t);
    elseif M >= 0.3 && M < 1
        Cl_a_w_inc = 2*pi*0.4/(1+2*0.7/A_w);
        Cl_a_t_inc = 2*pi*0.4/(1+2*0.7/A_t);
        Cl_a_w = Cl_a_w_inc/sqrt(1-M^2); % Correzione di Prandtl-Glauert per effetti di compressibilità
        Cl_a_t = Cl_a_t_inc/sqrt(1-M^2); % Correzione di Prandtl-Glauert per effetti di compressibilità
    elseif M == 1
         Cl_a_w_inc = 2*pi*0.50/(1+2*0.8/A_w);
         Cl_a_t_inc = 2*pi*0.50/(1+2*0.8/A_t);
         Cl_a_w = Cl_a_w_inc/sqrt(1-0.92^2); 
         Cl_a_t = Cl_a_t_inc/sqrt(1-0.92^2); 
    elseif M > 1 && M <= 1.1
        Cl_a_w = 4*0.65/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 1.1 && M <= 1.2
        Cl_a_w = 4*0.72/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 1.2 && M <= 1.3
        Cl_a_w = 4*0.74/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 1.3 && M <= 1.4
        Cl_a_w = 4*0.76/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 1.4 && M <= 1.5
        Cl_a_w = 4*0.78/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 1.5 && M <= 1.6
        Cl_a_w = 4*0.80/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 1.6 && M <= 1.7
        Cl_a_w = 4*0.81/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 1.7 && M <= 1.8
        Cl_a_w = 4*0.82/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 1.8 && M <= 1.9
        Cl_a_w = 4*0.83/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 1.9 && M <= 2
        Cl_a_w = 4*0.84/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2 && M <= 2.1
        Cl_a_w = 4*0.85/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2.1 && M <= 2.2
        Cl_a_w = 4*0.86/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2.2 && M <= 2.3
        Cl_a_w = 4*0.87/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2.3 && M <= 2.4
        Cl_a_w = 4*0.88/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2.4 && M <= 2.5
        Cl_a_w = 4*0.89/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2.5 && M <= 2.6
        Cl_a_w = 4*0.90/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2.6 && M <= 2.7
        Cl_a_w = 4*0.91/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2.7 && M <= 2.8
        Cl_a_w = 4*0.92/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2.8 && M <= 2.9
        Cl_a_w = 4*0.93/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 2.9 && M <= 3
        Cl_a_w = 4*0.94/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3 && M <= 3.1
        Cl_a_w = 4*0.95/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3.1 && M <= 3.2
        Cl_a_w = 4*0.96/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3.2 && M <= 3.3
        Cl_a_w = 4*0.97/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3.3 && M <= 3.4
        Cl_a_w = 4*0.98/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3.4 && M <= 3.5
        Cl_a_w = 4*0.99/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3.5 && M <= 3.6
        Cl_a_w = 4*1.00/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3.6 && M <= 3.7
        Cl_a_w = 4*1.01/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3.7 && M <= 3.8
        Cl_a_w = 4*1.02/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3.8 && M <= 3.9
        Cl_a_w = 4*1.03/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 3.9 && M <= 4
        Cl_a_w = 4*1.04/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4 && M <= 4.1
        Cl_a_w = 4*1.05/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4.1 && M <= 4.2
        Cl_a_w = 4*1.06/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4.2 && M <= 4.3
        Cl_a_w = 4*1.08/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4.3 && M <= 4.4
        Cl_a_w = 4*1.10/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4.4 && M <= 4.5
        Cl_a_w = 4*1.12/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4.5 && M <= 4.6
        Cl_a_w = 4*1.14/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4.6 && M <= 4.7
        Cl_a_w = 4*1.16/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4.7 && M <= 4.8
        Cl_a_w = 4*1.18/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4.8 && M <= 4.9
        Cl_a_w = 4*1.20/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 4.9 && M <= 5
        Cl_a_w = 4*1.22/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5 && M <= 5.1
        Cl_a_w = 4*1.24/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5.1 && M <= 5.2
        Cl_a_w = 4*1.26/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5.2 && M <= 5.3
        Cl_a_w = 4*1.28/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5.3 && M <= 5.4
        Cl_a_w = 4*1.30/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5.4 && M <= 5.5
        Cl_a_w = 4*1.32/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5.5 && M <= 5.6
        Cl_a_w = 4*1.34/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5.6 && M <= 5.7
        Cl_a_w = 4*1.36/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5.7 && M <= 5.8
        Cl_a_w = 4*1.38/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5.8 && M <= 5.9
        Cl_a_w = 4*1.40/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 5.9 && M <= 6
        Cl_a_w = 4*1.42/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6 && M <= 6.1
        Cl_a_w = 4*1.44/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6.1 && M <= 6.2
        Cl_a_w = 4*1.46/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6.2 && M <= 6.3
        Cl_a_w = 4*1.48/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6.3 && M <= 6.4
        Cl_a_w = 4*1.50/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6.4 && M <= 6.5
        Cl_a_w = 4*1.52/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6.5 && M <= 6.6
        Cl_a_w = 4*1.56/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6.6 && M <= 6.7
        Cl_a_w = 4*1.58/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6.7 && M <= 6.8
        Cl_a_w = 4*1.60/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6.8 && M <= 6.9
        Cl_a_w = 4*1.62/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 6.9 && M <= 7
        Cl_a_w = 4*1.64/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 7 && M <= 7.1
        Cl_a_w = 4*1.68/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 7.1 && M <= 7.2
        Cl_a_w = 4*1.70/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 7.2 && M <= 7.3
        Cl_a_w = 4*1.72/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 7.3 && M <= 7.4
        Cl_a_w = 4*1.74/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    elseif M > 7.4 
        Cl_a_w = 4*1.76/sqrt(M^2-1);
        Cl_a_t = Cl_a_w;
    end
end

end




