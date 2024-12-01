clear all
clc
close all

%% INPUT:
Mach_v = 0.3:0.1:8;
alpha_v = 0:5:50;
h = 11000;

% GEOMETRY:
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
S_wing = 3;        % m^2

GEO_funzione2 = struct('x', x, 'a', a, 'phi', phi, 'nose_type', nose_type, 'nose_data', data, 'S_wing', S_wing);
% x --> vettore delle posizioni longitudinali
% a --> verrore dei fìdiametri corrispettivi a tali posizioni
% phi --> angolo di rotazione attorno all'asse x
% nose_type --> 'C' per conical, 'TO' per tangent ogive
% nose_data --> database aerodinamica naso
% S_wing --> superficie portante totale (ala + tail)
% ,



%% TEST della funzione 1:
% tic
% [CL, CD] = file_funzione1(Mach_v, alpha_v, h, GEO_funzione1)
% toc
% figure
% hold on
% plot(Mach_v, CL, 'b')
% plot(Mach_v, CD, 'k')
% xlabel('Mach')
% ylabel('CD (nero), CL (blu)')


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
% save('data_CD_GEO1_S3m^2.mat', 'CD_mat')
% save('data_CL_GEO1_S3m^2.mat', 'CL_mat')


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


% Xcg = 6;
% tic
% [CD, CL_NKP, l_Cp_results, Cm] = file_funzione2_jeno(Mach_v, alpha_v, h, Xcg, GEO_funzione1)
% toc
% 
% % figure
% plot(Mach_v, CL_NKP, '.-', Color='b', LineWidth=2)
% title('CL vs mach')
% 
% % figure
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

Xcg = 6;

% Matrice di risultati:
h_v = 11000:500:84000;
CL_mat = zeros(length(Mach_v), length(alpha_v), length(h_v));       % righe: Mach da 0.3 a 8 passo 0.1
                                                                    % colonne:
                                                                    % angolo d'attacco da 0 a 50 passo 5
                                                                    % pagine: quote da 11000 a 84000 passo 500

CD_mat = zeros(length(Mach_v), length(alpha_v), length(h_v));
Xcp_mat = zeros(length(Mach_v), length(alpha_v), length(h_v));
tic
for p=1:length(h_v)

    h = h_v(p);
    [CD, CL, l_Cp_results, Cm] = file_funzione2_final(Mach_v, alpha_v, h, Xcg, GEO_funzione2)

    CL_mat(:, :, p) = CL;
    CD_mat(:, :, p) = CD;
    Xcp_mat(:, :, p) = l_Cp_results;
    p

end
toc

save('data_CD_GEO2_S3m^2.mat', 'CD_mat')
save('data_CL_GEO2_S3m^2.mat', 'CL_mat')
save('data_Xcp_GEO2_S3m^2.mat', 'Xcp_mat')

figure
plot(Mach_v, CL, '.-', Color='b', LineWidth=2)
title('CL vs mach')













