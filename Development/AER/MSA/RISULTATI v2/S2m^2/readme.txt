D_ogiva = 0.6 * 2;
D_base = 0.6 * 2;
L_ogiva = 2.4;
L_corpo_sup = 5;
L_spalla = 1;
L_corpo_inf = 12.85;
a = [0, D_ogiva, D_ogiva, D_base, D_base];
x = [0, L_ogiva, L_ogiva+L_corpo_sup, L_ogiva+L_corpo_sup+L_spalla, L_ogiva+L_corpo_sup+L_spalla+L_corpo_inf];
b = a;      % radius is constant (no ellipse)
% nel caso la sezione fosse un'ellisse:
phi = 0;       % IN CASE of ELLIPSE --> orientation wrt the normal velocity
data = xlsread('Dataset Cdn.xlsx', 'Default Dataset');
nose_type = 'C';    % C --> conical, TO --> tangent ogive
% INPUT Area dell'ala:
S_wing = 2;        % m^2
