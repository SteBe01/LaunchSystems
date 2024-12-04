function [CLAMP]=Handling(MASS,GEOMETRY,b)

% Teoria linea elastica appoggio-appoggio

% Ipotizzo lanciatore==trave e forze su baricentro
FS=1.5;
m = MASS.m_tot;
l=GEOMETRY.L_tot;
E=MASS.E;
J = MASS.J;

g0 = 9.80665; %m/s^2
F = m*g0;
a = l-b;

delta = @(x)  ( (-F*b*x^3) + F*b*(l^2 - b^2)*x )/(6*l*E*J);

phi_a = a*b*(l+b)*F/(6*E*l*J);
phi_b = a*b*(l+a)*F/(6*E*l*J);

CLAMP.phi_a = phi_a;
CLAMP.phi_b = phi_b;
x_vec = linspace(0,l,1000);

for i=1:length(x_vec)

f(i,1) = delta(x_vec(i));
end
CLAMP.f = f;
f_max = max(f);
CLAMP.f_max = f_max;

if f_max>l/3000

    fprintf('Too high displacement /n');

end
phi_max = max([abs(phi_a);abs(phi_b)]);
CLAMP.phi_max=phi_max;
if phi_max>(10^-3)

    fprintf('Too high rotations /n');

end



%% Forces on the clamp:

mu = 0.7; % Shigley's Mechanical Engineering Design: rubber-Al

CLAMP.MAT = material_selection(5); % Al 7XXX

R_a = FS*F*b/l;
R_b = FS*F*a/l;

F_clamp_a = R_a/mu;
F_clamp_b = R_b/mu;

F_clamp = max([F_clamp_a;F_clamp_b]);

A_clamp = F_clamp/CLAMP.MAT.sy;

CLAMP.A=A_clamp;

CLAMP.a = a;
CLAMP.b = b;
%R_LV =GEOMETRY.RLV;
%CLAMP.L_clamp = A_clamp/(2*pi*R_LV);








end