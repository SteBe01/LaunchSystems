function [Launcher_Prel_Des,Baseline] = Preliminary_Design(IS,EPSILON,Delta_V, M_pay,lambda0,Plot_switch)

options = optimoptions('fsolve', 'Display', 'none');
g0 = 9.80665; %m/s^2

TR.M_pay_max_mission = 400; % [kg] constarint on max volume

% Computation of error of our model wrt Baseline data:

P.M01 = 23130; % [kg]
[MP_model,~,~] = initialmass_opt(290.2*g0,289.4*g0,287.4*g0,0.083562229,0.095830454,0.140625,500,8165.975324,lambda0);

P.Err_M01_P = (abs(MP_model.M01 - P.M01))/P.M01; % error of our model wrt Real Pegasus data

L1.M01 = 25854; % [kg]
[ML1_model,~,~] = initialmass_opt_2stage(309*g0,328*g0,0.06,0.06,500,8462.129169,lambda0);

L1.Err_M01_L1 = (abs(ML1_model.M01 - L1.M01))/L1.M01; % error of our model wrt Real LauncherOne data

TR.Delta_V_tot = Delta_V;

TR.M_pay = M_pay; 

[TR.N_stages,~] = size(IS);

switch TR.N_stages

    case 2

c1 = IS(1,1)*g0;
c2 = IS(2,1)*g0;
eps_1 = EPSILON(1,1);
eps_2 = EPSILON(2,1);


[TR.M,TR.n,lambda_a] = initialmass_opt_2stage(c1,c2,eps_1,eps_2,TR.M_pay,TR.Delta_V_tot,lambda0);


 case 3

c1 = IS(1,1)*g0;
c2 = IS(2,1)*g0;
c3 = IS(3,1)*g0;
eps_1 = EPSILON(1,1);
eps_2 = EPSILON(2,1);
eps_3 = EPSILON(3,1);

[TR.M,TR.n,lambda_a] = initialmass_opt(c1,c2,c3,eps_1,eps_2,eps_3,TR.M_pay,TR.Delta_V_tot,lambda0);

end



%% Fairing Volume Evaluation:

% Empirically compute the volume of the fairing by taking data from the
% datasheets and approximating the volume to a cylinder + truncated cone +
% small cylinder for the nose: 
% Use base diameter, nose diameter and the heights of the 3 different parts

% Inizialization of all the data from datasheets (state page and where are
% found) and then use the function V_fair_approx to compute the
% approximated volume and empty one

% Electron: payload user guide 7.0, pg 18,23; Carbon composite fairing

El.D_fair_base_int_El = 1.070; % [m]
El.D_fair_base_ext_El = 1.200; % [m]
El.H_fair_base_El = 0.566 + 0.385 + 0.0705; % [m]
El.D_fair_conetrap_int_El = 0.2783;
El.D_fair_conetrap_ext_El = 0.2783 + (0.13);
El.H_fair_conetrap_int_El = 1.3586; 
El.H_fair_nose_El = 0.119;
El.M_fair_El = 44; % [kg]
El.M_pay_max_El = 300;% [kg]
El.fair_thick =  (1.200-1.070)/2;

El.V_fair_El = V_fair_approx(El.D_fair_base_int_El,El.D_fair_base_ext_El,El.H_fair_base_El,El.D_fair_conetrap_int_El,El.D_fair_conetrap_ext_El,El.H_fair_conetrap_int_El,El.H_fair_nose_El);

% LauncherOne: service guide august 2020, pg 14, 15; Carbon composite fairing

L1.D_fair_base_int_L1 = 1.2624; % [m]
L1.D_fair_base_ext_L1 = 1.2624 + (2*0.11176); % [m]
L1.H_fair_base_L1 = 2.1234 + (4*0.00635); % [m]
L1.D_fair_conetrap_int_L1 = 0.44196;
L1.D_fair_conetrap_ext_L1 = 0.44196 + (2*0.11176);
L1.H_fair_conetrap_int_L1 = 3.5433-2.1234; 
L1.H_fair_nose_L1 = 3.63 - 3.5433;
L1.M_fair_L1 = 145;
L1.M_pay_max_L1 = 500; % [kg]
L1.fair_thick = 0.11176;

L1.V_fair_L1 = V_fair_approx(L1.D_fair_base_int_L1,L1.D_fair_base_ext_L1,L1.H_fair_base_L1,L1.D_fair_conetrap_int_L1,L1.D_fair_conetrap_ext_L1,L1.H_fair_conetrap_int_L1,L1.H_fair_nose_L1);

% Pegasus: Pegasus user guide 1; pg 13,15,40  
% aluminum honeycomb core with graphite/epoxy skins for the cylindrical and ogive sections of the fairing and a monocoque graphite/epoxy nose cap

P.D_fair_base_int_P = 1.153; % [m]
P.D_fair_base_ext_P = 1.27; % [m]
P.H_fair_base_P = 1.11 + 0.4; % [m]
P.D_fair_conetrap_int_P = 0.709;
P.D_fair_conetrap_ext_P = (1.27-1.153) + 0.709;
P.H_fair_conetrap_int_P = 2.656 - 1.11 + 0.4; 
P.H_fair_nose_P = 3.63 - 3.5433;
P.M_fair_P = 170; % [kg] Maggi's ex
P.M_pay_max_P = 443; % [kg]
P.fair_thick = (1.27-1.153)/2;

P.V_fair_P = V_fair_approx(P.D_fair_base_int_P,P.D_fair_base_ext_P,P.H_fair_base_P,P.D_fair_conetrap_int_P,P.D_fair_conetrap_ext_P,P.H_fair_conetrap_int_P,P.H_fair_nose_P);


El.rho_pay_El = El.M_pay_max_El/El.V_fair_El.V_fair_tot;
L1.rho_pay_L1 = L1.M_pay_max_L1/L1.V_fair_L1.V_fair_tot;
P.rho_pay_P = P.M_pay_max_P/P.V_fair_P.V_fair_tot;

El.rho_pay_El_empty = El.M_pay_max_El/El.V_fair_El.V_fair_tot_int;
L1.rho_pay_L1_empty = L1.M_pay_max_L1/L1.V_fair_L1.V_fair_tot_int;
P.rho_pay_P_empty = P.M_pay_max_P/P.V_fair_P.V_fair_tot_int;

El.rho_fair_El = El.M_fair_El/(El.V_fair_El.V_fair_tot - El.V_fair_El.V_fair_tot_int);
L1.rho_fair_L1 = L1.M_fair_L1/(L1.V_fair_L1.V_fair_tot - L1.V_fair_L1.V_fair_tot_int);
P.rho_fair_P = P.M_fair_P/(P.V_fair_P.V_fair_tot - P.V_fair_P.V_fair_tot_int);

% Volume and mass of fairing of our launcher:

% Compute volume and mass of fairing of our launcher via interpolation of
% data


% Interpolation of data for whole volume of fairing:
x_rhop_eq = [L1.M_pay_max_L1,P.M_pay_max_P,El.M_pay_max_El];
y_rhop_eq = [L1.rho_pay_L1,P.rho_pay_P,El.rho_pay_El];
V_rhop_eq = polyfit(x_rhop_eq,y_rhop_eq,1);
y_rhop_line = @(x) V_rhop_eq(1)*x + V_rhop_eq(2);

% Interpolation of data for empty volume of fairing:
x_rhope_eq = [L1.M_pay_max_L1,P.M_pay_max_P,El.M_pay_max_El];
y_rhope_eq = [L1.rho_pay_L1_empty,P.rho_pay_P_empty,El.rho_pay_El_empty];
V_rhope_eq = polyfit(x_rhope_eq,y_rhope_eq,1);
y_rhope_line = @(x) V_rhope_eq(1)*x + V_rhope_eq(2);

% Use of data interpolation to get values of our fairing
V_fair = TR.M_pay_max_mission/(V_rhop_eq(1)*TR.M_pay_max_mission + V_rhop_eq(2));
TR.Fair.V_fair = V_fair;
V_fair_empty = TR.M_pay_max_mission/(V_rhope_eq(1)*TR.M_pay_max_mission + V_rhope_eq(2)); % [m^3] Empty volume
TR.Fair.V_fair_empty = V_fair_empty;
% Mass of fairing computed from mean densities of our baselines:

M_fair = mean([L1.rho_fair_L1;El.rho_fair_El;P.rho_fair_P]) * (V_fair - V_fair_empty); % [kg]
TR.Fair.M_fair = M_fair;
%% Volume of launcher: From Excel on Baseline: 

% Compute volume and empty one of our launcher using baseline data, the
% necessary data i all initialised in terms of baseline and the values from
% TR.M are used which were taken from the optimization_stage_2stages

% Definition of values for baseline:

L1.fn_ratio1 = 12; % [-] fineness ratio of LauncherOne
L1.fn_ratio2 = 14; % [-] fineness ratio of LauncherOne
L1.fn_avg = mean([L1.fn_ratio1;L1.fn_ratio2]);
L1.V_empty = (0.356*(9.31+1.87)*pi*(((1.8+1.5)/2)^2)/4); % [m^3] Maggi ex
L1.rho_empty = L1.M01/L1.V_empty; % [kg/m^3]
L1.V_fair_L1.Ltot = 3.63; % [m]
L1.fn_nose = L1.V_fair_L1.Ltot/L1.D_fair_base_ext_L1; % [-] Fineness ratio of nose of LauncherOne
L1.rho_launcher = 614.9798737; % [kg/m^3]
L1.V_launcher = 42.04040019; % [m^3]

P.fn_ratio = 13.31; % [-] fineness ratio of Pegasus
P.V_empty = 0.805*pi*(1.27^2)/4 + (1.905*pi*(1.27^2)/4); % [m^3]
P.rho_empty = P.M01/P.V_empty; % [kg/m^3]
P.V_fair_P.Ltot = 2.65; % [m]
P.fn_nose = P.V_fair_P.Ltot/P.D_fair_base_ext_P; % [-] Fineness ratio of nose of Pegasus
P.rho_launcher = 1041.647642; % [kg/m^3]
P.V_launcher = 21.40839099; % [m^3]

El.M01 = 12550;
El.fn_ratio = 14;
El.V_launcher = 19.226547; % [m^3]
El.rho_launcher = El.M01/El.V_launcher;
El.V_empty = 4.08099; % [m^3]
El.rho_empty = El.M01/El.V_empty;
El.V_fair_El.Ltot = 2.5; % [m]
El.fn_nose = El.V_fair_El.Ltot/El.D_fair_base_ext_El; % [-]

% Interpolate values from the 3 baselines as before:

% Interpolation of full densities (fairing does not count as it is hyp full)
y_rhol_eq = [L1.rho_launcher;P.rho_launcher;El.rho_launcher];
x_rhol_eq = [L1.M01;P.M01;El.M01];
V_rhol_eq = polyfit(x_rhol_eq,y_rhol_eq,1);
y_rhol_line = @(x) V_rhol_eq(1)*x + V_rhol_eq(2);

% Volume of launcher is obtained from interpolation
TR.V_launcher = TR.M.M01/(V_rhol_eq(1)*TR.M.M01 + V_rhol_eq(2));

% Interpolation of empty densities (fairing does not count as it is hyp full)
y_rhole_eq = [L1.rho_empty;P.rho_empty;El.rho_empty];
x_rhole_eq = [L1.M01;P.M01;El.M01];
V_rhole_eq = polyfit(x_rhole_eq,y_rhole_eq,1);
y_rhole_line = @(x) V_rhole_eq(1)*x + V_rhole_eq(2);

% Empty volume of launcher is obtained from interpolation
TR.V_launcher_empty = TR.M.M01/(V_rhole_eq(1)*TR.M.M01 + V_rhole_eq(2));

% TR.V_launcher = V_massl_eq(1)*TR.M.M01 + V_massl_eq(2);
% TR.V_launcher_empty = V_massle_eq(1)*TR.M.M01 + V_massle_eq(2);

% Effective volume as difference between whole and empty volume
TR.V_launcher_effective = TR.V_launcher - TR.V_launcher_empty; % [m^3]

% Interpolation of fineness ratios:
y_fn_eq = [mean([L1.fn_ratio1;L1.fn_ratio2]);P.fn_ratio;El.fn_ratio];
x_fn_eq = [L1.M01;P.M01;El.M01;];
V_fn_eq = polyfit(x_fn_eq,y_fn_eq,1);

%TR.fn_ratio = P.fn_ratio; % [-] L/D=f (Pegasus Baseline) imposed ??
TR.fn_ratio = round(V_fn_eq(1)*TR.M.M01 + V_fn_eq(2));

TR.Diameter = ((4*TR.V_launcher)/(pi*TR.fn_ratio) )^(1/3);  % [m]
TR.Length = TR.fn_ratio * TR.Diameter; % [m] Whole body length

fair_thick = mean([El.fair_thick;P.fair_thick;L1.fair_thick]);
TR.Fair.fair_thick = fair_thick;
% TR.fair_length = TR.V_fair/( ((pi/4)*(TR.Diameter^2)) + ((pi/3)*(((TR.Diameter^2)/4)+(((TR.Diameter-2*TR.fair_thick)^2)/4)+((TR.Diameter*(TR.Diameter-2*TR.fair_thick))/4))) + ((pi/4)*((TR.Diameter-2*TR.fair_thick)^2)));
fair_length = (V_fair*4)/(pi*(TR.Diameter^2)); % [m] Fairing length
TR.Fair.fair_length = fair_length;
TR.body_length = TR.Length - fair_length; % [m] Body length

Baseline.El = El;
Baseline.L1 = L1;
Baseline.P = P;

Launcher_Prel_Des = TR;


if Plot_switch == 1

switch TR.N_stages

    case 2

f = @(x) c1*log(c1*x -1) + c2*log(c2*x -1) - log(x)*(c1+c2) - c1*log(c1*eps_1) - c2*log(c2*eps_2) - TR.Delta_V_tot;
lambda_vec = ((10^-4):0.000001:(10^-3));
lambda_vec = lambda_vec';
y = zeros(length(lambda_vec),1);
for i = 1:length(lambda_vec)
y(i) = abs(f(lambda_vec(i)));
end

figure()
plot(lambda_vec,y);
hold on
plot(lambda_a,f(lambda_a),'o');
xlabel('\lambda',Interpreter='tex');
ylabel('Cost Function');
legend('Cost Function','Fsolve solution');

  case 3
f = @(x) c1*log(c1*x -1) + c2*log(c2*x -1) + c3*log(c3*x -1) - log(x)*(c1+c2+c3) - c1*log(c1*eps_1) - c2*log(c2*eps_2) - c3*log(c3*eps_3) - TR.Delta_V_tot;

lambda_vec = ((10^-4):0.000001:(10^-3));
lambda_vec = lambda_vec';
y = zeros(length(lambda_vec),1);
for i = 1:length(lambda_vec)
y(i) = abs(f(lambda_vec(i)));
end

figure()
plot(lambda_vec,y);
hold on
plot(lambda_a,f(lambda_a),'o');
xlabel('\lambda',Interpreter='tex');
ylabel('Cost Function');
legend('Cost Function','Fsolve solution');



end
    

% Body plot:

figure() 
plot([0 fair_length TR.Length],[0 TR.Diameter TR.Diameter]/2,'k',[0 fair_length TR.Length],-[0 TR.Diameter TR.Diameter]/2,'k');
axis equal
grid on
title('Geometry');
xlabel('Length [m]');
ylabel('Radius [m]');

% Plot to see where we are wrt Baseline

x_vec = linspace(0,1000,300);
figure()
plot(x_vec,y_rhop_line(x_vec));
hold on;
plot(x_rhop_eq,y_rhop_eq,'o');
hold on;
plot(TR.M_pay_max_mission,(V_rhop_eq(1)*TR.M_pay_max_mission + V_rhop_eq(2)),'^');
xlabel('Maximum Payload Mass $[kg]$',Interpreter='latex');
ylabel('Density of Fairing $[\frac{kg}{m^3}]$',Interpreter='latex');
title('Linear Interpolation of Payload Mass and Fairing Density');
legend('','Baseline','Team rocket');

x_vec = linspace(0,1000,300);
figure()
plot(x_vec,y_rhope_line(x_vec));
hold on;
plot(x_rhope_eq,y_rhope_eq,'o');
hold on;
plot(TR.M_pay_max_mission,(V_rhope_eq(1)*TR.M_pay_max_mission + V_rhope_eq(2)),'^');
xlabel('Maximum Payload Mass $[kg]$',Interpreter='latex');
ylabel('Density of Empty Fairing $[\frac{kg}{m^3}]$',Interpreter='latex');
title('Linear Interpolation of Payload Mass and Fairing Empty Density');
legend('','Baseline','Team rocket');

x_vec = linspace(0,50000,300);
figure()
plot(x_vec,y_rhol_line(x_vec));
hold on;
plot(x_rhol_eq,y_rhol_eq,'o');
hold on;
plot(TR.M.M01,(V_rhol_eq(1)*TR.M.M01 + V_rhol_eq(2)),'^');
xlabel('GLOM $[kg]$',Interpreter='latex');
ylabel('Density of Launcher $[\frac{kg}{m^3}]$',Interpreter='latex');
title('Linear Interpolation of GLOM and Launcher Density');
legend('','Baseline','Team rocket');

x_vec = linspace(0,50000,300);
figure()
plot(x_vec,y_rhole_line(x_vec));
hold on;
plot(x_rhole_eq,y_rhole_eq,'o');
hold on;
plot(TR.M.M01,(V_rhole_eq(1)*TR.M.M01 + V_rhole_eq(2)),'^');
xlabel('GLOM $[kg]$',Interpreter='latex');
ylabel('Density of Empty Launcher $[\frac{kg}{m^3}]$',Interpreter='latex');
title('Linear Interpolation of GLOM and Launcher Empty Density');
legend('','Baseline','Team rocket');


end


end