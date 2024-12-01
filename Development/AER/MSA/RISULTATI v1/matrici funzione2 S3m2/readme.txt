
h_v = 11000:500:84000;
Mach_v = 0.3:0.1:8;
alpha_v = 0:5:50;


% GEOMETRY:
D_ogiva = 0.525 * 2;
D_base = 0.7 * 2;
L_ogiva = 2.1;
L_corpo_sup = 4.292;
L_spalla = 1.751;
L_corpo_inf = 12.013;
a = [0, D_ogiva, D_ogiva, D_base, D_base];
x = [0, L_ogiva, L_ogiva+L_corpo_sup, L_ogiva+L_corpo_sup+L_spalla, L_ogiva+L_corpo_sup+L_spalla+L_corpo_inf];
b = a;      % radius is constant (no ellipse)
% nel caso la sezione fosse un'ellisse:
phi = 0;       % IN CASE of ELLIPSE --> orientation wrt the normal velocity
data = xlsread('Dataset Cdn.xlsx', 'Default Dataset');
nose_type = 'C';    % C --> conical, TO --> tangent ogive
% INPUT Area dell'ala:
S_wing = 3;        % m^2
