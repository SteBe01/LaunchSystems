% CALCOLA ALCUNE COSE AL VARIARE DELLA MANETTA, COMPRESA TRA IL 50 E IL 100%
% E SALVA I DATI IN "data_throttle.mat" PER POTERLI UTILIZZARE ALTROVE
% SENZA DOVER CHIAMARE OGNI VOLTA IL CEAM

clear
clc
close all
addpath("CEA");
options  = optimoptions("fsolve", "MaxFunctionEvaluations",1000, 'Display','none');

g0       = 9.81;
OF       = 2.6;
epsilonc = 3.3611;
epsilon  = 330.5028;
A_t      = 0.0028;

throttle          = struct();
throttle.pp_ch    = linspace(55/2, 55, 1000);
throttle.TT       = NaN(1,length(throttle.pp_ch));
throttle.m_dot    = NaN(1,length(throttle.pp_ch));
throttle.p_e      = NaN(1,length(throttle.pp_ch));
throttle.m_dot_fu = NaN(1,length(throttle.pp_ch));
throttle.m_dot_ox = NaN(1,length(throttle.pp_ch));

for i = 1:length(throttle.TT)
    data     = fromCEA(OF, epsilonc, epsilon, throttle.pp_ch(i), 0);
    throttle.p_e(i)      = data.p_e;
    throttle.m_dot(i)    = throttle.pp_ch(i)*1e5*A_t/data.c_star;
    throttle.m_dot_fu(i) = throttle.m_dot(i)*(1/(1+OF));
    throttle.m_dot_ox(i) = throttle.m_dot(i)*(OF/(1+OF));
    throttle.TT(i)       = data.i_sp*g0*throttle.m_dot(i);
end
throttle.manetta = throttle.TT/throttle.TT(end)*100;

save("data_throttle.mat", "throttle")

%%

function out = fromCEA(OF, epsilonc, epsilon, Pc, flag)
    output = CEA('problem', 'rocket', 'equilibrium','fr','nfz',2,'o/f',OF,'subsonic(ae/at)',epsilonc,'supsonic(ae/at)',epsilon, ...
    'p(bar)',Pc,'reactants', ...
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