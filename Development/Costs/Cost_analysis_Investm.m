%% Cost Analysis
% Fonte Manned Spacecraft, Forza
clc; clear
% Investment 
n_dev = 4;        % 4 anni di development
D = 10000*n_dev;  % Non-recurring development cost
Stage1 = 2000*1400;
B747 = 6.5*10^6;
Motori1 = 8*450000;
Ltot = 10;
CFUr = Motori1 + B747 + Stage1;
Invs = CFUr + D;
Price = 4.5*10^6;
Soldi0 = -Invs;
for i = 1:Ltot

L = i;            % number of launches
N = 1;            % years for the amortiz. (lifetime of the project)
% Cost development coefficients (Liquid Rocket)
n = L ;           % total number of units
S = 0.7;            % learning constant

CFUsp = 2000*200 +450000 ;

B = 1 + log(S)/log(2);   % learning factor

Stage2= ((CFUsp)* n^B)/L;

Cops = 10^5;

F1 = 0.15;        % 0.10-0.15 aging factor
F2 = 1.10;        % 1.05-1.10 aging factor

if L >=2
Cf = F1*F2*(Stage1+Motori1)/Ltot;
else
    Cf = 0;
end
Cr = 30000;

%Ci = 0.15*Cp;     % insurance, "the insurance rate has typically been considered 
                  % as a percentage, say around 15%, of the cost of the launch

Cll = (Cr + Cf + Stage2 + Cops)/(1-0.15);


Soldi0 = Soldi0 - Cll + Price;
Soldi(i) = Soldi0;


end

plot([1:Ltot],Soldi/10^6)
hold on; xlabel('Number of launches'); ylabel('Mln $'); grid on
plot([1:Ltot],Soldi/10^6,'o')
title('ROI, 4.5 Mln per launch')
yline(0)


               



