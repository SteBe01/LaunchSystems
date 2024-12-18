clear all
clc
close all

%% INPUT:
Mach_v = 3:1:8;
alpha_v = 0:5:10;           % deg
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
nose_type = 'TO';    % C --> conical, TO --> tangent ogive
% INPUT Area dell'ala
S_wing = 0;        % m^2

GEO_funzione2 = struct('x', x, 'a', a, 'phi', phi, 'nose_type', nose_type, 'nose_data', data, 'S_wing', S_wing);
% x --> vettore delle posizioni longitudinali
% a --> verrore dei fìdiametri corrispettivi a tali posizioni
% phi --> angolo di rotazione attorno all'asse x
% nose_type --> 'C' per conical, 'TO' per tangent ogive
% nose_data --> database aerodinamica naso
% S_wing --> superficie portante totale (ala + tail)
% ,



%% TEST della funzione 2:

Xcg = 6;
tic
[CD, CL_NKP, l_Cp_results, Cm] = file_funzione2_final(Mach_v, alpha_v, h, Xcg, GEO_funzione2)
toc

% figure
% plot(Mach_v, CL_NKP, '.-', Color='b', LineWidth=2)
% title('CL vs mach')
% 
% figure
% plot(Mach_v, CD, '.-', Color='k', LineWidth=2)
% title('CD vs mach')
% 
% figure
% plot(Mach_v, Cm, 'k')
% title('Cm vs Mach')
% 
% figure
% plot(Mach_v, l_Cp_results, 'b')
% title('Xcp vs Mach')
% 
% Xcg = 6;
% 
% % Matrice di risultati:
% h_v = h;
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



%% TEST pannelli:

mainPanel



%% PLOT:

figure
hold on
grid minor
plot(Mach_v, CL_NKP, '.-', LineWidth=2)
plot(Mach_v, CL_ModNew, '--', LineWidth=2)
title('CL su Mach')
legend('Jorgensen alpha 0°', 'Jorgensen alpha 5°', 'Jorgensen alpha 10°', 'Modified Newton 0°', 'Modified Newton 5°', 'Modified Newton 10°')
ylim([0, 2])

figure
hold on
grid minor
plot(Mach_v, CD*0.8, '.-', LineWidth=2)
plot(Mach_v, CD_ModNew*1.7, '--', LineWidth=2)
legend(...
    '$0^\circ$ Semi-Empirical', ...
    '$5^\circ$ Semi-Empirical', ...
    '$10^\circ$ Semi-Empirical', ...
    '$0^\circ$ Modified Newton', ...
    '$5^\circ$ Modified Newton', ...
    '$10^\circ$ Modified Newton', ...
    'Interpreter', 'latex', 'FontSize', 12);
ylim([0, 1])
dim = [0.2 0.5 0.3 0.1]; % [x, y, larghezza, altezza] in unità normalizzate (0-1)
    annotation('textbox', dim, 'String', 'ISA altitude = 20000 m', ...
           'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'Interpreter','latex');
xlabel('Mach', 'Interpreter','latex', LineWidth=1.5)
ylabel('Cd [-]', 'Interpreter','latex', LineWidth=1.5)
ax = gca; % Ottieni l'oggetto Axes corrente
ax.TickLabelInterpreter = 'latex';







