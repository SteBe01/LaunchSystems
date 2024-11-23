clear
close all
clc
addpath("CEA");

%% INPUT

%DATA
gimbal = 5; %degs
D_th   = 0.06;  % 7 cm extern
D_ch   = 0.11; % 12.6 cm ext
D_ex_old = 0.24; % 25 cm ext
OF        = 2.6;
g0=9.81;

%VARIABLES
d_rocket=1.20;
N_engines=7;
L_engine=0.89;
h_launch = 11900;

%% PRELIMINARY COMPUTING
epsilonc = (D_ch/D_th)^2;
P_amb = p_amb_fun(h_launch);
A_th = pi*D_th^2/4;

%% RUTHERFORD VALIDATION

P_ch=validateEngine(D_th, D_ch, D_ex_old);

%% NEW CONFIGURATION MAXIMIZING EPSILON

epsilon = configureEngines (d_rocket, N_engines, L_engine, gimbal, D_th, 1);
Ae = epsilon * A_th;

% can safely handle Pe/Pamb of 0.35, verify we are above
CEAout = fromCEA(OF, epsilonc, epsilon, P_ch, 0);
Pe = CEAout.p_e;
cstar = CEAout.c_star;
Isp = CEAout.i_sp;

Pressure_ratio = Pe/P_amb;
if Pressure_ratio<0.35
    error("Too much overexpansion")
end

mdot_launch = P_ch*A_th/cstar;
T_single_engine = mdot_launch*Isp*g0 + (Pe-P_amb)*Ae;
T_liftoff = N_engines*T_single_engine;

z_opt_km = z_of_pressure(Pe)/1e3;

zz=11900:100:30000;

pp=p_amb_fun(zz);
TT=N_engines*mdot_launch*Isp*g0 + (Pe-pp)*Ae;
figure
plot(zz/1e3,TT/1e3)
ylabel("Thrust [KN]")
xlabel("Altitude [km]")
yyaxis right
plot(zz/1e3,pp/1e3)
ylabel("Pressure [kPa]")




%% HELPERS

function P_ch = validateEngine(D_th, D_ch, D_ex)

A_th       = pi*D_th^2/4;
OF        = 2.6;
epsilonc  = (D_ch/D_th)^2;
epsilon   = (D_ex/D_th)^2;
T_sl      = 24e3;

P_ch = fzero(@(Pc) MassFlowDifference(OF, epsilonc, epsilon, Pc, T_sl, A_th),5.7357e6);

end

%%
function err = MassFlowDifference(OF, epsilonc, epsilon, Pc, T, At)
CEAout = fromCEA(OF, epsilonc, epsilon, Pc, 0);

P_amb = 101325;
g0 = 9.81;
Ae = epsilon*At;

Isp = CEAout.i_sp;
Pe  = CEAout.p_e;
cstar = CEAout.c_star;

mdot_launch = (T - (Pe-P_amb)*Ae)/(Isp*g0);
mdot_cea = Pc*At/cstar;
err=mdot_cea-mdot_launch;
end

%%
function out = fromCEA(OF, epsilonc, epsilon, Pc, flag)
output = CEA('problem', 'rocket', 'equilibrium','fr','nfz',2,'o/f',OF,'subsonic(ae/at)',epsilonc,'supsonic(ae/at)',epsilon, ...
    'p(bar)',Pc*1e-5,'reactants', ...
    'fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100.,'t(k)',298.15,...
    'oxid','O2(L)','O',2,'wt%',100, 't(k)',90.17,'output','transport','mks','end');

if flag
    out = output.output.froz.pressure(end) * 1e5;
else
    out = struct('p_e', output.output.froz.pressure(end)*1e5,...
        'i_sp', output.output.froz.isp(end),...
        'i_sp_void', output.output.froz.isp_vac(end),...
        'c_star', output.output.froz.cstar(end),...
        'gamma', output.output.froz.gamma(end),...
        'temp',  output.output.froz.temperature);
end
end

%%
function p = p_amb_fun(z)


p0 = 101325;
p = p0*(1-2.25577e-5*z).^5.25588;

end

function z = z_of_pressure(p)
    p0 = 101325;
    z = -log(p/p0)*8600;
end