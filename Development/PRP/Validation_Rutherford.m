clear;clc;

A_t       = 3.2e-3;
OF        = 2.6;
Pc        = 5.8e6;
epsilonc  = 3.36;
epsilon   = 16;
A_e       = epsilon*A_t;

output = CEA('problem', 'rocket', 'equilibrium','fr','nfz',2,'o/f',OF,'subsonic(ae/at)',epsilonc,'supsonic(ae/at)',epsilon, ...
    'p(bar)',Pc*1e-5,'reactants',  ...
    'fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100.,'t(k)',298.15,...
    'oxid','O2(L)','O',2,'wt%',100, 't(k)',90.17,'output','transport','mks','end');