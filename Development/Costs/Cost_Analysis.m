%% Cost Analysis
% Fonte Manned Spacecraft, Forza
clc; clear
L_tot = 10;       % Total number of launcher to be considered
Price = 4.5;      % [mln dollars]
for i = 1:L_tot
n_dev = 4;        % 4 anni di development
L = i;            % number of launches
D = 10000*n_dev;  % Non-recurring development cost
N = 1;            % years for the amortiz. (lifetime of the project)
% Cost development coefficients (Liquid Rocket)
n = L ;           % total number of units
S = 0.7;            % learning constant
ai = 187;
bi = 0.60; 

Cd = (D)/L ;      % senza ammortamento, per aggiungerlo si complica la formula

%CFUr = 77*10^6;   % Reused theoretical cost of the first unit 
Stage1 = 2000*1400;
B747 = 6.5*10^6;
Motori1 = 8*450000;

CFUr = Motori1 + B747 + Stage1;
%CFUsp = 50*10^6;  % theoretical cost of the first unit of second stage 
CFUsp = 2000*200 +450000 ;

B = 1 + log(S)/log(2);   % learning factor

fl = n^(B-1);     % learning curve factor
Costsecondstage= ((CFUsp)* n^B)/L/CFUsp;
P = CFUr  + (CFUsp)* n^B; % total production of n vehicles

Cp = P/L;  

%Cops = 10^7;      % Average for unmanned and reusable
Cops = 10^5;

F1 = 0.15;        % 0.10-0.15 aging factor
F2 = 1.10;        % 1.05-1.10 aging factor

if i >= 2
Cf = F1*F2*(Stage1+Motori1)/L_tot;
else
    Cf = 0;
end

Cr = 30000;

Ci = 0.15*Cp;     % insurance, "he insurance rate has typically been considered 
                  % as a percentage, say around 15%, of the cost of the launch

Cll = Ci + Cr +Cf + Cp + Cops;
Cl(i) = Cll/10^6;
end

plot(1:L_tot,Cl)
hold on
plot(1:L_tot,Cl,'o')
yline(Price,'--','Price we want to sell a launch')
[~,ROI_ID] = max(Cl<4.5);
grid on
xline(ROI_ID,'--','Min N. of launches to break even')
hold on; xlabel('Number of launches'); ylabel('Mln $'); grid on
title('Cost per launch to have ROI')



% 360 kg + 40 kg (40 a motore), secondo stadio, massa totale inerte 1400
% 450 mila a motore

               



