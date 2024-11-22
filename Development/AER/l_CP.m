function [l_CP] = l_CP(M, alpha, l_s,V_S,r_N,l_w,r, s_w,cr_w, q_inf, m_w, d, K_BW, K_WB, K_N_W, LW, s_t, l_t,...
    cr_t,m_t, K_BT, K_TB, K_N_T, LT,Cl_N,Cl_WB, Cl_BW, Cl_BT, Cl_TB, Cl_TV)

% INPUTS:
% M = mach number
% alpha = angle of attack
% l_s = distance from most forward point of body to shoulder of body  nose [m]
% V_S = volume of body nose up to should [m^3]
% r_N = body radius at shoulder of nose [m]
% l_w = distance from most forward point of body to intersection of wing
%       leading edge and body [m]
% r = body radiu [m]
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


% CENTRE OF PRESSURE OF BODY NOSE

l_N = l_s*(1-V_S/(pi*r_N^2*l_s));



% CENTRE OF PRESSURE OF WING IN PRESENCE OF BODY

x_cr_WB_a = 1/(1-r/s_w)*(2*(1/3+r^4/s_w^4)*atan(s_w/r)+2/3*r^3/s_w^3*log(((s_w^2+r^2)/(2*s_w^2))^2*s_w/r)-1/3*r^3/s_w^3*...
    (2*pi-1+s_w^2/r^2))/((1+r^2/s_w^2)^2*atan(s_w/r)-r^2/s_w^2*(pi+(s_w/r-r/s_w)))-r/s_w/(1-r/s_w);

l_WB = l_w + cr_w*x_cr_WB_a;



% CENTRE OF PRESSURE ON BODY DUE TO WING

beta = sqrt(abs(M^2-1));

if M > 1

    M_BW = 4*q_inf*alpha*m_w/(3*pi*beta)*cr_w^3*(sqrt(1+2*beta*d/cr_w)*((2*m_w*beta+5)/(3*(m_w*beta+1)^2)+beta*d/cr_w/(3*...
        (m_w*beta+1))-(beta*d/cr_w)^2/(beta*m_w))+1/sqrt(m_w^2*beta^2-1)*((1+beta*d/cr_w)^3-(beta*d/cr_w)^3/(m_w^2*beta^2)-1/...
        (1+m_w*beta)^2)*acos((1+beta*d/cr_w*(m_w*beta+1))/(m_w*beta+beta*d/cr_w*(m_w*beta+1)))+(beta*d/cr_w)^3*1/(m_w^2*beta^2)*...
        acosh(1+cr_w/(beta*d))-((2*m_w*beta+5)/(3*(m_w*beta+1)^2))-(1-(1/(m_w*beta+1))^2)/sqrt(m_w^2*beta^2-1)*acos(1/(m_w*beta)));

elseif M < 1

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

l_TB = l_t + cr_t*x_cr_TB_a;



% CENTRE OF PRESSURE ON BODY DUE TO TAIL


if M > 1

    M_BT = 4*q_inf*alpha*m_t/(3*pi*beta)*cr_t^3*(sqrt(1+2*beta*d/cr_t)*((2*m_t*beta+5)/(3*(m_t*beta+1)^2)+beta*d/cr_t/(3*...
        (m_t*beta+1))-(beta*d/cr_t)^2/(beta*m_t))+1/sqrt(m_t^2*beta^2-1)*((1+beta*d/cr_t)^3-(beta*d/cr_t)^3/(m_t^2*beta^2)-1/...
        (1+m_t*beta)^2)*acos((1+beta*d/cr_t*(m_t*beta+1))/(m_t*beta+beta*d/cr_t*(m_t*beta+1)))+(beta*d/cr_t)^3*1/(m_t^2*beta^2)*...
        acosh(1+cr_t/(beta*d))-((2*m_t*beta+5)/(3*(m_t*beta+1)^2))-(1-(1/(m_t*beta+1))^2)/sqrt(m_t^2*beta^2-1)*acos(1/(m_t*beta)));

elseif M < 1

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

l_CP = (l_N*Cl_N++l_WB*Cl_WB+l_BW*Cl_BW+l_BT*Cl_BT+l_TB*Cl_TB+l_TV*Cl_TV)/(Cl_N+Cl_WB+Cl_BW+Cl_BT+Cl_TB+Cl_TV);














