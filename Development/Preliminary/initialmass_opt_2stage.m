function [M,n,lambda] = initialmass_opt_2stage(c1,c2,eps_1,eps_2,M_pay,Delta_V,lambda0)

options = optimoptions('fsolve', 'Display', 'none');


f = @(x) c1*log(c1*x -1) + c2*log(c2*x -1) - log(x)*(c1+c2) - c1*log(c1*eps_1) - c2*log(c2*eps_2) - Delta_V;

lambda = fsolve(f,lambda0,options);

n.n1 = (c1*lambda - 1)/(c1*eps_1*lambda); % [-] 1/MR1 % must be larger than 1
n.n2 = (c2*lambda - 1)/(c2*eps_2*lambda); % [-] 1/MR2

M.M2 = ((n.n2 -1)/(1 - n.n2*eps_2)) * (M_pay); % [kg] Mass of stage 2
M.M1 = ((n.n1 -1)/(1 - n.n1*eps_1)) * (M.M2+ M_pay); % [kg] Mass of stage 1

M.Ms1 = eps_1*M.M1; % [kg] Mass of structure of 1
M.Ms2 = eps_2*M.M2; % [kg] Mass of structure of 2

M.Mp1 = M.M1 - M.Ms1; % [kg] Mass of propellaft of 1
M.Mp2 = M.M2 - M.Ms2; % [kg] Mass of propellant of 2

M.M02 = M_pay+ M.M2; % [kg] Mass of stack 2
M.M01 = M.M02 + M.M1; % [kg] Mass of stack 1


M.MR1 = 1/n.n1;
M.MR2 = 1/n.n2;

if n.n1<1 | isreal(n.n1)==0
    M.M01 = NaN;
     n.n1 = NaN;
     M.MR1 = NaN;
elseif n.n2<1 | isreal(n.n2)==0
 M.M01 = NaN;
  n.n2 = NaN;
  M.MR2 = NaN;

end

end