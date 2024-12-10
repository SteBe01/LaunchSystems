
%% run this with the selected architecture (r = row of the matrix)
close all;


r = 34;

Is = [324.5; 343.0]; %[s] 
g0 = 9.8065; %[m/s^2]
m_pay_eff = 251; %[kg] effective payload mass
m_adapter = 0.0755 * m_pay_eff + 50; %[kg] adapter mass
m_pay = m_pay_eff + m_adapter; %[kg] mass of payload+adapter

%recover reference data:
mpay0 = M_it(r).pay_effective; %[kg]
dv0 = M_it(r).dv; %[km/s]

%get masses of the structures
str1 = M_it(r).str1; 
str2 = M_it(r).str2; 

%get masses of the propellants
prp1 = M_it(r).M0 - M_it(r).M0end;
prp2 = M_it(r).M1 - M_it(r).M1end;

%get m0, m0end, m1, m1end
m0 = str1 + str2 + prp1 + prp2 + m_pay; %[kg]
m0end = m0 - prp1; %[kg]
m1 = m0end - str1; %[kg]
m1end = m1 - prp2; %[kg]

%get mr1, mr2
mr1 = m0/m0end; %[-]
mr2 = m1/m1end; %[-]

%dv estimation
dv1 = Is(1)*g0 * log(mr1); %[m/s]
dv2 = Is(2)*g0 * log(mr2); %[m/s]
dv = (dv1 + dv2)/1e3; %[km/s]

%tradeoff between dv and m_pay:
ddv_dmpay = (dv-dv0)/(m_pay-mpay0); %[km/(s*kg)]






