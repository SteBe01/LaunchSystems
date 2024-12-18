clear all
clc
close all

%% INPUT:
Mach_v = 0.5:0.05:8;
alpha_v = 0;           % deg
alpha_v = deg2rad(alpha_v);     % rad
h = 20000;

% GEOMETRY:
D_ogiva = 1.2;
D_base = 1.2;
L_ogiva = 2.4;
L_corpo_sup = 4.4040;
L_spalla = 0;
L_corpo_inf = 19.5503 - L_spalla - L_corpo_sup - L_ogiva;
a = [0, D_ogiva, D_ogiva, D_base, D_base];
x = [0, L_ogiva, L_ogiva+L_corpo_sup, L_ogiva+L_corpo_sup+L_spalla, L_ogiva+L_corpo_sup+L_spalla+L_corpo_inf];
b = a;      % radius is constant (no ellipse)
% nel caso la sezione fosse un'ellisse:
phi = 0;       % IN CASE of ELLIPSE --> orientation wrt the normal velocity
data = xlsread('Dataset Cdn.xlsx', 'Default Dataset');
nose_type = 'C';    % C --> conical, TO --> tangent ogive
% INPUT Area dell'ala
S_wing = 1;        % m^2

GEO_funzione2 = struct('x', x, 'a', a, 'phi', phi, 'nose_type', nose_type, 'nose_data', data, 'S_wing', S_wing);
% x --> vettore delle posizioni longitudinali
% a --> verrore dei fìdiametri corrispettivi a tali posizioni
% phi --> angolo di rotazione attorno all'asse x
% nose_type --> 'C' per conical, 'TO' per tangent ogive
% nose_data --> database aerodinamica naso
% S_wing --> superficie portante totale (ala + tail)
% ,



%% TEST della funzione 1:
tic
[CL1, CD1] = file_funzione1(Mach_v, alpha_v, h, GEO_funzione2)
toc
figure
hold on
plot(Mach_v, CL1, 'b')
plot(Mach_v, CD1, 'k')
xlabel('Mach')
ylabel('CD (nero), CL (blu)')


% % Matrice di risultati:
% h_v = 11000:500:84000;
% CL_mat = zeros(length(Mach_v), length(alpha_v), length(h_v));       % righe: Mach da 0.3 a 8 passo 0.1
%                                                                     % colonne:
%                                                                     % angolo d'attacco da 0 a 50 passo 5
%                                                                     % pagine: quote da 11000 a 84000 passo 500
% 
% CD_mat = zeros(length(Mach_v), length(alpha_v), length(h_v));
% tic
% for p=1:length(h_v)
% 
%     h = h_v(p);
%     [CL, CD] = file_funzione1(Mach_v, alpha_v, h, GEO_funzione1);
% 
%     CL_mat(:, :, p) = CL;
%     CD_mat(:, :, p) = CD;
%     p
% 
% end
% toc
% 
% save('data_CD_GEO1_S1m^2.mat', 'CD_mat')
% save('data_CL_GEO1_S1m^2.mat', 'CL_mat')


%% TEST della funzione 2:

% Mach_v = 0.1:0.1:8;
% alpha_v = 0:10;
% h = 15000;
% 
% % GEOMETRY:
% D_ogiva = 0.525 * 2;
% D_base = 0.7 * 2;
% L_ogiva = 2.1;
% L_corpo_sup = 4.292;
% L_spalla = 1.751;
% L_corpo_inf = 12.013;
% a = [0, D_ogiva, D_ogiva, D_base, D_base];
% x = [0, L_ogiva, L_ogiva+L_corpo_sup, L_ogiva+L_corpo_sup+L_spalla, L_ogiva+L_corpo_sup+L_spalla+L_corpo_inf];
% b = a;      % radius is constant (no ellipse)
% % nel caso la sezione fosse un'ellisse:
% phi = 0;       % IN CASE of ELLIPSE --> orientation wrt the normal velocity
% data = xlsread('Dataset Cdn.xlsx', 'Default Dataset');
% nose_type = 'C';    % C --> conical, TO --> tangent ogive
% % INPUT Area dell'ala:
% S_wing = 0;        % m^2
% 
% GEO_funzione1 = struct('x', x, 'a', a, 'phi', phi, 'nose_type', nose_type, 'nose_data', data, 'S_wing', S_wing)


Xcg = 6;
tic
[CD, CL_NKP, l_Cp_results, Cm] = file_funzione2_final(Mach_v, alpha_v, h, Xcg, GEO_funzione2)
toc

% figure
% hold on
% plot(Mach_v, CL_NKP, '.-', Color='b', LineWidth=2)
% plot(Mach_v, CL1, '--', Color='k', LineWidth=2)
% title('CL vs mach')
% 
% figure
% hold on
% plot(Mach_v, CD, '.-', Color='b', LineWidth=2)
% plot(Mach_v, CD1, '--', Color='k', LineWidth=2)
% title('CD vs mach')
% 
% figure
% plot(Mach_v, Cm, 'k')
% title('Cm vs Mach')
% 
figure
hold on
grid minor
plot(Mach_v, l_Cp_results, color = "#0072BD", LineWidth=2)
xlabel('Mach [-]', 'FontName', 'Palatino Linotype', 'FontSize', 12, LineWidth=2)
ylabel('X_c_p [m]', 'FontName', 'Palatino Linotype', 'FontSize', 12, LineWidth=2)
ax = gca; % Ottieni l'oggetto Axes corrente
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 12;



%% REPORT:
figure
hold on
grid minor
plot(rad2deg(alpha_v), CL_NKP, LineWidth=1.5)
xlabel('Alpha [°]')
ylabel('Cl [-]')
legend('Mach 0.5', 'Mach 1.0', 'Mach 1.5', 'Mach 2.0', 'Mach 2.5', 'Mach 3.0', 'Mach 3.5', 'Mach 4.0')
% title('Lift coefficent Cl [-]')

figure
hold on
grid minor
plot(Mach_v, CD*0.85, LineWidth=2.25)
xlabel('Mach')
ylabel('Cd [-]')
legend('AoA 0°', 'AoA 2°', 'AoA 4°', 'AoA 6°', 'AoA 8°', 'AoA 10°')
% title('Drag Coefficient [-]')





% Xcg = 6;
% 
% % Matrice di risultati:
% h_v = 11000:500:84000;
% CL_mat = zeros(length(Mach_v), length(alpha_v), length(h_v));       % righe: Mach da 0.3 a 8 passo 0.1
%                                                                     % colonne:
%                                                                     % angolo d'attacco da 0 a 50 passo 5
%                                                                     % pagine: quote da 11000 a 84000 passo 500
% 
% CD_mat = zeros(length(Mach_v), length(alpha_v), length(h_v));
% Xcp_mat = zeros(length(Mach_v), length(alpha_v), length(h_v));
% tic
% for p=1:length(h_v)
% 
%     h = h_v(p);
%     [CD, CL, l_Cp_results, Cm] = file_funzione2_final(Mach_v, alpha_v, h, Xcg, GEO_funzione2);
% 
%     CL_mat(:, :, p) = CL;
%     CD_mat(:, :, p) = CD;
%     Xcp_mat(:, :, p) = l_Cp_results;
%     p
% 
% end
% toc
% 
% save('data_CD_GEO3_S1m^2.mat', 'CD_mat')
% save('data_CL_GEO3_S1m^2.mat', 'CL_mat')
% save('data_Xcp_GEO3_S1m^2.mat', 'Xcp_mat')
% 
% figure
% hold on
% plot(Mach_v, CL, '.-', Color='b', LineWidth=2)
% plot(Mach_v, CD, '.-', Color='k', LineWidth=2)
% title('CL vs mach')














