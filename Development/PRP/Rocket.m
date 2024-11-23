clear; clc;
addpath("CEA");

g0 = 9.80665;

d_t       = 0.06;
A_t       = pi*d_t^2/4;
OF        = 2.6;
p_chamber = 5.8e6;
epsilonc  = 3.4;

h0       = 13000;
options  = optimoptions("fsolve", "MaxFunctionEvaluations",1000, 'Display','none');
deltap   = @(epsilon) fromCEA(OF, epsilonc, epsilon, p_chamber*1e-5, 1)-p_amb(h0);
epsilon  = fsolve(deltap, 31, options);
out_prop = fromCEA(OF, epsilonc, epsilon, p_chamber*1e-5, 0);
A_e      = epsilon*A_t;
d_e      = sqrt(4*A_e/pi);
A_c      = A_t*epsilonc;
d_c      = sqrt(4*A_c/pi);
L_car    = 0.62;
l_ch     = L_car*A_t/A_c;

% T_req = ;   % thrust profile
% h     = ;   % altitude profile
% 
% deltaT = @(p_ch) thrust(h, OF, epsilonc, epsilon, p_ch, A, m_dot) - T_req;
% p_c = fsolve(deltaT, p_chamber, options);


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

function p = p_amb(z)
    p0 = 101325;
    p = p0*exp(-z/8600);
end

function T = thrust(h, OF, epsilonc, epsilon, Pc, A, m_dot)

    output = CEA('problem', 'rocket', 'equilibrium','fr','nfz',2,'o/f',OF,'subsonic(ae/at)',epsilonc,'supsonic(ae/at)',epsilon, ...
    'p(bar)',Pc,'reactants', ...
    'fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100.,'t(k)',298.15,...
    'oxid','O2(L)','O',2,'wt%',100, 't(k)',90.17,'output','transport','mks','end');

    T_static = A*(output.output.froz.pressure(end)-p_amb(h));
    T = output.output.froz.isp(end)*g0*m_dot + T_static;

end