clear, clc
close all
clear dyn
[Out] = Trajectory();

Max_Q = max(Out.qdyn);
q = Max_Q;
t_max_Q = find(Max_Q==Out.qdyn);
nx = Out.acc(t_max_Q,1);
nx = nx(1,1);
T = Out.stages.stg1.Thrust;
Cd = Out.stages.stg1.Cd;
A_ref = pi*(Out.stages.stg1.d^2/4);
m = Out.mass(t_max_Q);
g = Out.gravity(t_max_Q);

D = Cd*q*A_ref;

P_balance = T - D - m*g*(nx + 1);

if P_balance> 1

nx = (T-D)/(m_tot*g0) - 1;

end



Max_M = max(Out.moment);

t_max_M = find(Max_M==Out.moment);

% figure()
% plot(T,moment)
% hold on;
% plot(T,qdyn)
% legend('Moment','Max_Q');