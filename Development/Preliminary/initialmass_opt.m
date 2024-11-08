function [M,n,lambda] = initialmass_opt(c1,c2,c3,eps_1,eps_2,eps_3,M_pay,Delta_V,initial_guess)
options = optimoptions('fsolve', 'Display', 'none');

f = @(x) c1*log(c1*x -1) + c2*log(c2*x -1) + c3*log(c3*x -1) - log(x)*(c1+c2+c3) - c1*log(c1*eps_1) - c2*log(c2*eps_2) - c3*log(c3*eps_3) - Delta_V;

lambda = fsolve(f,initial_guess,options);
%lambda = fzero(f,4.6e-04); % lambda vicino 0.1 e check per local min

n.n1 = (c1*lambda - 1)/(c1*eps_1*lambda); % [-] 1/MR1 % must be larger than 1
n.n2 = (c2*lambda - 1)/(c2*eps_2*lambda); % [-] 1/MR2
n.n3 = (c3*lambda - 1)/(c3*eps_3*lambda); % [-] 1/MR3

M.M3 = ((n.n3 -1)/(1 - n.n3*eps_3)) * (M_pay); % [kg] Mass of stage 3
M.M2 = ((n.n2 -1)/(1 - n.n2*eps_2)) * (M.M3 + M_pay); % [kg] Mass of stage 2
M.M1 = ((n.n1 -1)/(1 - n.n1*eps_1)) * (M.M2 + M.M3 + M_pay); % [kg] Mass of stage 1

M.Ms1 = eps_1*M.M1; % [kg] Mass of structure of 1
M.Ms2 = eps_2*M.M2; % [kg] Mass of structure of 2
M.Ms3 = eps_3*M.M3; % [kg] Mass of structure of 3

M.Mp1 = M.M1 - M.Ms1; % [kg] Mass of propellaft of 1
M.Mp2 = M.M2 - M.Ms2; % [kg] Mass of propellant of 2
M.Mp3 = M.M3 - M.Ms3; % [kg] Mass of propellant of 3

M.M03 = M.M3 + M_pay; % [kg] Mass of stack 3
M.M02 = M.M03 + M.M2; % [kg] Mass of stack 2
M.M01 = M.M02 + M.M1; % [kg] Mass of stack 1

M.MR1 = 1/n.n1;
M.MR2 = 1/n.n2;
M.MR3 = 1/n.n3;

if n.n1<1 | isreal(n.n1)==0
    M.M01 = NaN;
     n.n1 = NaN;
     M.MR1 = NaN;
elseif n.n2<1 | isreal(n.n2)==0
 M.M01 = NaN;
  n.n2 = NaN;
  M.MR2 = NaN;
elseif n.n3<1 | isreal(n.n3)==0
     M.M01 = NaN;
     n.n3 = NaN;
     M.MR3 = NaN;

end

end