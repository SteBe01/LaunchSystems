%%%%%%%%%%%%%%% INPUT Geometria %%%%%%%%%%%%%%%%%%%

% INPUT diametri e posizioni corrispettive
x = [0 3 8 9 15];
a = [0 1.5 1.5 1.8 1.8];
b = a;      % radius is constant (no ellipse)

% nel caso la sezione fosse un'ellisse:
phi = 0;       % IN CASE of ELLIPSE --> orientation wrt the normal velocity
a_max = max(a);
b_max = max(b);


% INPUT forma dell'ogiva
ln = x(2);
dn = a(2);
nose_type = 'C';    % C --> conical, TO --> tangent ogive
% Nose data loading from excel file:
data = excel_load;

% INPUT Area dell'ala:
A_w = 1;        % m^2

%% REF. DATA:

% HP.: sezione circolare
A_b = pi * a_max^2;     % area della base
A_r = A_b;              % area di riferimento
A_p = trapz(x, a)*2;    % area planform



%% Parametri per NKP e l_Cp
r_N = 0.6; 
S_w = (1.92+0.764)*1.08;
s_w = 1.08+0.65;
r = 0.65;
betaA_w = 0.8*(2*s_w)^2/S_w;
lambda_w = 0.5;
m_w = cot(deg2rad(15));
cr_w = 1.92;
r_w = 0.65;
c_w = 1.3;
S_t = (1.92+0.764)*1.08;
s_t = 1.08+r;
betaA_t = 0.8*(2*s_t)^2/S_t;
lambda_t = 0.5;
m_t = cot(deg2rad(15));
cr_t = 1.92;
A_w = (2*s_w)^2/S_w;
A_t = (2*s_t)^2/S_t;
r_t = 0.65;
l_w = 6.5;
l_t = 15-cr_t;
c_t = 1.3;
l_s = 3;
V_S = pi*r_N^2*l_s/3;
d = 1.3;

S = A_p + S_w + S_t;
S_surface = S_w + S_t;



%% PLOT GEOMETRIA:
% figure
% hold on
% grid minor
% plot(x, a/2, 'k', x, -a/2,'k', LineWidth=1.75)
% plot([x(end), x(end), x(end)-2*A_w/sqrt(3*A_w)], [a(end)/2, a(end)/2+sqrt(3*A_w), a(end)/2], 'b', LineWidth=2)
% plot([x(end), x(end), x(end)-2*A_w/sqrt(3*A_w)], [-a(end)/2, -a(end)/2-sqrt(3*A_w), -a(end)/2], 'b', LineWidth=2)
% plot([x(end), x(end)-2*A_w/sqrt(4*A_w)], [0, 0], 'b', LineWidth=2)
% axis equal
% grid on
% set(gca, 'FontSize', 30)
% xlabel('x [m]', FontSize=35)
% ylabel('y [m]', FontSize=35)
% title('Launcher shape', FontSize=40)


%% Excel file loading:
function data = excel_load

% Excel file loading:
data = xlsread('Dataset Cdn.xlsx', 'Default Dataset');

end
