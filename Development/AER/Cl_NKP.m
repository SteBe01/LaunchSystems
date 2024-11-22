function [Cl_NKP] = Cl_NKP(alpha, M, h, r_N, S_w, s_w, r, betaA_w, lambda_w, m_w, cr_w, r_w, c_w, Cl_w,S_t,s_t,betaA_t,...
                           lambda_t,m_t,cr_t,A_t,V_inf,r_t, l_w,l_t, x_hat_t,c_t, Cl_a_w, Cl_a_t)

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



d = 2*r;  % d = body diameter [m]

if h > 84000
    Cl_NKP = 0;
end

alpha = deg2rad(alpha); % [rad]

% deal with AoA
if alpha <= deg2rad(90) && alpha >= 0
    alpha = + alpha;
elseif alpha <= deg2rad(180) && alpha >= deg2rad(90)
    alpha = deg2rad(180) - alpha;
end



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

beta = sqrt(abs(M^2-1));

if betaA_w*(1+lambda_w)*(1/(beta*m_w)+1) < 4 % slender-body theory
    if r/s_w == 0 % the combination is all wing
        K_BW = 0;
    elseif r/s_w >= 1 % there is a very small exposed wing
        K_BW = 2;
    else
        K_BW = (1-r^2/s_w^2)^2 - 2/pi * ((1+r^4/s_w^4)*(1/2*atan(1/2*(s_w/r-r/s_w))+pi/4)-r^2/s_w^2*((s_w/r-r/s_w)+2*atan(r/s_w)))/(1-r/s_w)^2;
    end
elseif betaA_w*(1+lambda_w)*(1/(beta*m_w)+1) >= 4
    if M > 1
        K_BW = 8*beta*m_w/(pi*sqrt(beta^2*m_w^2-1)*(1+lambda_w)*(beta*d/cr_w)*(s_w/r-1)*beta*Cl_a_w)*((beta*m_w/(1+beta*m_w))*...
            ((beta*m_w+1)*beta*d/cr_w+beta*m_w/(beta*m_w))^2*acos(1+(1+beta*m_w)*beta*d/cr_w/(beta*m_w+(beta*m_w+1)*beta*d/cr_w))+...
            sqrt(beta^2*m_w^2-1)/(beta*m_w+1)*(sqrt(1+s_w*beta*d/cr_w)-1)-sqrt(beta^2*m_w^2-1)/(beta*m_w)*(beta*d/cr_w)^2*...
            acosh(1+cr_w/(beta*d))-beta*m_w/(1+beta*m_w)*acos(1/(beta*m_w)));
    elseif M < 1
        K_BW = 16*(beta*m_w/(1+m_w*beta))^2/(pi*(1+lambda_w)*(beta*d/cr_w)*(s_w/r-1)*beta*Cl_a_w)*((beta*m_w+(1+m_w*beta)*...
            beta*d/cr_w/(beta*m_w))^(3/2)+(beta*m_w+(1+m_w*beta)*beta*d/cr_w/(beta*m_w))^(1/2)-2-((1+m_w*beta)*beta*d/cr_w/...
            (m_w*beta))^2*atanh(sqrt(beta*m_w/(beta*m_w+(1+m_w*beta)*beta*d/cr_w))));
    end
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

if betaA_t*(1+lambda_t)*(1/(beta*m_t)+1) < 4 % slender-body theory
    if r/s_t == 0 % the combination is all tail
        K_BT = 0;
    elseif r/s_t >= 1 % there is a very small exposed tail
        K_BT = 2;
    else
        K_BT = (1-r^2/s_t^2)^2 - 2/pi * ((1+r^4/s_t^4)*(1/2*atan(1/2*(s_t/r-r/s_t))+pi/4)-r^2/s_t^2*((s_t/r-r/s_t)+2*atan(r/s_t)))/(1-r/s_t)^2;
    end
elseif betaA_t*(1+lambda_t)*(1/(beta*m_t)+1) >= 4
    if M > 1
        K_BT = 8*beta*m_t/(pi*sqrt(beta^2*m_t^2-1)*(1+lambda_t)*(beta*d/cr_t)*(s_t/r-1)*beta*Cl_a_t)*((beta*m_t/(1+beta*m_t))*...
            ((beta*m_t+1)*beta*d/cr_t+beta*m_t/(beta*m_t))^2*acos(1+(1+beta*m_t)*beta*d/cr_t/(beta*m_t+(beta*m_t+1)*beta*d/cr_t))+...
            sqrt(beta^2*m_t^2-1)/(beta*m_t+1)*(sqrt(1+s_t*beta*d/cr_t)-1)-sqrt(beta^2*m_t^2-1)/(beta*m_t)*(beta*d/cr_t)^2*...
            acosh(1+cr_t/(beta*d))-beta*m_t/(1+beta*m_t)*acos(1/(beta*m_t)));
    elseif M < 1
        K_BT = 16*(beta*m_t/(1+m_t*beta))^2/(pi*(1+lambda_t)*(beta*d/cr_t)*(s_t/r-1)*beta*Cl_a_t)*((beta*m_t+(1+m_t*beta)*...
            beta*d/cr_t/(beta*m_t))^(3/2)+(beta*m_t+(1+m_t*beta)*beta*d/cr_t/(beta*m_t))^(1/2)-2-((1+m_t*beta)*beta*d/cr_t/...
            (m_t*beta))^2*atanh(sqrt(beta*m_t/(beta*m_t+(1+m_t*beta)*beta*d/cr_t))));
    end
end

Cl_BT = K_BT * Cl_a_t * alpha * S_t/S_w;



% LIFT ON TAIL SECTIONS DUE TO TAIL VORTICES

cl = 0.1; % giustificare con i paper
f_w = Cl_w*S_w/(2*cl*c_w);
f_t = Cl_t*S_t/(2*cl*c_t);
Gamma_m = V_inf*K_WB*alpha*Cl_a_w*S_w/(4*(f_w-r_w));
h_t = (l_t+x_hat_t-l_w-cr_w)*sin(alpha);
f_i = f_w*r_t^2/(f_w^2+h_t^2);
h_i = h_t*r_t^2/(f_w^2+h_t^2);
[L1] = L(lambda_t, r_t, s_t, f_w, h_t);
[L2] = L(lambda_t, r_t, s_t, -f_w, h_t);
[L3] = L(lambda_t, r_t, s_t, f_i, h_i);
[L4] = L(lambda_t, r_t, s_t, -f_i, h_i);

i = 2/(1+lambda_t)*(L1-L2-L3+L4);


Cl_TV = Cl_a_w*Cl_a_t*K_WB*alpha*i*(s_t-r_t)/(2*pi*A_t*(f_w-r_w));


% LIFT ON WING AFTERBODY DUE TO WING VORTICES

Cl_BV = -4*Gamma_m/(S_w*V_inf)*((f_w^2-r_w^2)/f_w-f_t+r_t^2/sqrt(f_t^2+r_t^2));



% TOTAL LIFT COEFFICIENT

Cl_NKP = Cl_N+Cl_WB+Cl_BW+Cl_TB+Cl_BT+Cl_TV+Cl_BV;


end















