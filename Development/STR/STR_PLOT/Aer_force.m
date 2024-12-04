function [AER]=Aer_force(q,alpha)
%% INPUT STR:

% INPUT diametri e posizioni corrispettive:
D_ogiva = 0.525 * 2;
D_base = 0.7 * 2;
L_ogiva = 2.1;
L_corpo_sup = 4.292;
L_spalla = 1.751;
L_corpo_inf = 12.013;

% Se lo voglio unito metto stesso diametro ogiva e base

% Input flusso:
%alpha = deg2rad(20);        % angolo d'attacco
M = 1.5;                    % mach
h = 11000;                  % altitude
[T,a_sound,P,rho,nu] = atmosisa(h);
%q = 0.5 * rho * (a_sound * M)^2;        % dynamic preassure

% INPUT Area ali:
A_t = 1;
A_w = 0;        % m^2
a = [0, D_ogiva, D_ogiva, D_base, D_base];
x = [0, L_ogiva, L_ogiva+L_corpo_sup, L_ogiva+L_corpo_sup+L_spalla, L_ogiva+L_corpo_sup+L_spalla+L_corpo_inf];
% 12.013 13.764 18.056 20.156]
% x = [0 3 8 10 18];
% a = [0 1.5 1.5 1.8 1.8];
b = a;      % radius is constant (no ellipse)

f_b = x(end) / max(a);
% % nel caso la sezione fosse un'ellisse:
% phi = 0;       % IN CASE of ELLIPSE --> orientation wrt the normal velocity
a_max = max(a);
b_max = max(b);

A_nose = pi * (a(2)^2)/4;


% REF. DATA:
% HP.: sezione circolare
A_b = pi * a_max^2;     % area della base
A_ref = A_b;              % area di riferimento
A_p = trapz(x, a)*2;    % area planform



%% Normal Forces:

% Cylinder1:
% CN_cyl1 = sin(2*alpha) * cos(alpha/2) + 1.3 * f_b * (sin(alpha))^2;
CN_cyl1 = 0;

% Cylinder2:
% CN_cyl2 = sin(2*alpha) * cos(alpha/2) + 1.3 * f_b * (sin(alpha))^2;
CN_cyl2 = 0;

% Nose:
CN_n_alpha = 2; %/ A_ref * ( A_nose );
CN_nose = CN_n_alpha * alpha;

% Shoulder:
% CN_sh_alpha = 2 / A_ref * ( A_sh_f - A_sh_i );
% CN_sh = CN_sh_alpha * alpha;
A_sh_f = pi * a(4)^2/4;
A_sh_i = pi * a(3)^2/4;
CN_sh_alpha = 8 / A_ref * ( A_sh_f - A_sh_i );
CN_sh = CN_sh_alpha * alpha;

% Tail:
if M^2 >= 1+(8/(pi*f_b))^2
   Cn_s = abs((4*abs(sin(alpha)*cos(alpha))/(M^2-1)^(1/2)+ 2*sin(alpha)^2)*(A_t/A_ref));
elseif M^2 < 1+(8/(pi*f_b))^2
   Cn_s = abs(((pi*f_b/2)*abs(sin(alpha)*cos(alpha))+2*sin(alpha)^2)*(A_t/A_ref));
end
CN_tail = Cn_s;

% Wing:
if M^2 >= 1+(8/(pi*f_b))^2
   Cn_s = abs((4*abs(sin(alpha)*cos(alpha))/(M^2-1)^(1/2)+ 2*sin(alpha)^2)*(A_w/A_ref));
elseif M^2 < 1+(8/(pi*f_b))^2
   Cn_s = abs(((pi*f_b/2)*abs(sin(alpha)*cos(alpha))+2*sin(alpha)^2)*(A_w/A_ref));
end
CN_wing = Cn_s;



%% Axial forces:

% SKIN FRICTION:
%%%%%%%%% % VAN DRIEST II: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if M == 0
    error('Mach = 0, non è possibile calcolare Cd_friction');
end
gamma = 1.4;
[T,a_sound,P,rho,nu] = atmosisa(h);
V = M * a_sound;                   % freestream velocity (m/s)

% Reynolds:
Re_x = V * max(a) / nu;        % Reynolds number based on length

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


% Nose:
f_nose = x(2) / a(2);
A_wet_nose = sqrt( (a(2)/2)^2 + (x(2))^2 ) * pi * a(2)/2;
Cf_nose = Cf * A_wet_nose/A_ref;
CA_fr_nose = Cf_nose * (1 + 0.5 / (f_nose) * A_wet_nose) / A_ref;

% Cyl1:
f_cyl1 = ( x(3)-x(2) ) / ( a(3) );
A_wet_cyl1 = 2 * a(3) * pi * ( x(3) - x(2) );
CA_fr_cyl1 =  A_wet_cyl1 * Cf / A_ref * (1 + 0.5 / (f_cyl1) * A_wet_cyl1) / A_ref;

% Shoulder:
f_sh = ( x(4)-x(3) ) * 2/ ( a(4) + a(3) );
A_wet_sh = sqrt( (a(4)/2)^2 + (x(4)-x(3)+x(2))^2 ) * pi * a(4)/2 - sqrt( (a(2)/2)^2 + (x(2))^2 ) * pi * a(2)/2;
CA_fr_sh =  A_wet_sh * Cf / A_ref * (1 + 0.5 / (f_sh) * A_wet_sh) / A_ref;

% Cyl2:
f_cyl2 = ( x(5)-x(4) ) / ( a(5) );
A_wet_cyl2 = 2 * a(5) * pi * ( x(5) - x(4) );
CA_fr_cyl2 =  A_wet_cyl2 * Cf / A_ref * (1 + 0.5 / (f_cyl2) * A_wet_cyl2) / A_ref;

% Tail:
A_wet_tail = A_t * 4;
CA_fr_tail = A_wet_tail * Cf / A_ref;

% Wing:
A_wet_wing = A_w * 2;
CA_fr_wing = A_wet_wing * Cf / A_ref;

%%%
CA_tail = CA_fr_tail;
CA_wing = CA_fr_wing;
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l = x(end);
d = max(a);
y = a;

fn_nose = x(2)/a(2);        % finesse ratio of nose cone

% Deal with AoA:
if alpha <= deg2rad(90) && alpha >= 0
    alpha = + alpha;
elseif alpha >= deg2rad(180) && alpha <= deg2rad(360)
    alpha = deg2rad(180);
end

%%%%%% WAVE: %%%%%%%%%%%%%
% Nose:
if alpha <= deg2rad(90)
    theta = atan(0.5 / fn_nose);
elseif alpha > deg2rad(90)
    theta = pi/2;
end
if M <= 1
    Ca_w_nose = 0.8 * sin(theta)^2;
elseif M > 1
    beta = sqrt(M^2-1);
    Ca_w_nose = (4 * sin(theta)^2 * (2.5 + 8 * beta * sin(theta))) / (1 + 16 * beta * sin(theta));
end
%%%
CA_nose = Ca_w_nose + CA_fr_nose;
%%%

% Shoulder:
if alpha <= deg2rad(90)
    theta = atan(0.5 / f_sh);
elseif alpha > deg2rad(90)
    theta = pi/2;
end
if M <= 1
    Ca_w = 0.8 * sin(theta)^2;
elseif M > 1
    beta = sqrt(M^2-1);
    Ca_w = (4 * sin(theta)^2 * (2.5 + 8 * beta * sin(theta))) / (1 + 16 * beta * sin(theta));
end
%%%
CA_sh = Ca_w - Ca_w_nose + CA_fr_sh;
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%% BASE PREASSURE: %%%%%%%%%%%%%%%%%%%%%
gamma = 1.4;
if M >= 1
    Cp_b = 2/(gamma * M^2) * ( (2/(gamma+1))^1.4 * (1/M)^2.8 * (2 * gamma * M^2 - gamma + 1)/(gamma + 1) - 1 );
    Ca_b = -Cp_b;
elseif M < 1
    Ca_b = 0.12 + 0.13 * M^2;
end
%%%
CA_cyl2 = Ca_b + CA_fr_cyl2;
%%%
%%%%%%%%%%%%%%%%%

%%%
CA_cyl1 = CA_fr_cyl1;
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% SUM of contributions:

% Nose:
cF_nose = [CA_nose, CN_nose] * 10;
pos_nose = [(x(2)+x(1))/2, 0];
F_nose = q * cF_nose;
AER.F_nose = norm(F_nose);
%fprintf('\nF_nose: %fN lungo x, %fN lungo y', F_nose(1), F_nose(2))

% Cyl1:
cF_cyl1 = [CA_cyl1, CN_cyl1] * 10;
pos_cyl1 = [(x(3)+x(2))/2, 0];
F_cyl1 = q * cF_cyl1;
%fprintf('\nF_upperbody: %fN lungo x, %fN lungo y', F_cyl1(1), F_cyl1(1))

% Shoulder:
cF_sh = [CA_sh, CN_sh] * 10;
pos_sh = [(x(4)+x(3))/2, 0];
F_sh = q * cF_sh;
%fprintf('\nF_shoulder: %fN lungo x, %fN lungo y', F_sh(1), F_sh(2))

% Cyl2:
cF_cyl2 = [CA_cyl2, CN_cyl2] * 10;
pos_cyl2 = [(x(5)+x(4))/2, 0];
F_cyl2 = q * cF_cyl2;
AER.D = norm(F_cyl2)+norm(F_cyl1);
%fprintf('\nF_lowerbody: %fN lungo x, %fN lungo y', F_cyl2(1), F_cyl2(2))

% Tail:
cF_tail = [CA_tail, CN_tail] * 10;
x_tail = x(end);
pos_tail = [x_tail, 0];
F_tail = q * cF_tail;
%fprintf('\nF_tail: %fN lungo x, %fN lungo y', F_tail(1), F_tail(2))

% Wing:
cF_wing = [CA_wing, CN_wing] * 10;
x_wing = x(end)/2;
pos_wing = [x_wing, 0];
F_wing = q * cF_wing;
AER.F_tail = norm(F_tail)+norm(F_wing);
%fprintf('\nF_wing: %fN lungo x, %fN lungo y', F_wing(1),  F_wing(2))

% Rename variables:
F_upperbody = F_cyl1;
F_shoulder = F_sh;
F_lowerbody = F_cyl2;





end