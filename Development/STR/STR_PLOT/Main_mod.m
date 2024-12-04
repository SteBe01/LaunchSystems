% TODO
% Implement the fuel derivative as done during classes

%% Simulator

clear, clc
close all
clear dyn

[stages, params, init] = loadMission();

% Compute masses and t_burn total
fn = fieldnames(stages);
for ii = 1:length(fn)
    stages.(fn{ii}).m_prop = stages.(fn{ii}).m0 * (1 - 1/stages.(fn{ii}).MR);
    stages.(fn{ii}).m_dot = stages.(fn{ii}).Thrust / (stages.(fn{ii}).Isp * params.g0);
    stages.(fn{ii}).t_burn_tot = stages.(fn{ii}).m_prop / stages.(fn{ii}).m_dot;
    if ii < numel(fn)
        stages.(fn{ii+1}).m0 = stages.(fn{ii}).m0 - stages.(fn{ii}).m_prop;
    end
end

% Trajectory propagation
y0_stg1 = [init.x0 init.z0 init.vx0 init.vz0 init.theta0 init.thetaDot0];

options_stg1 = odeset('RelTol',1e-8, 'MaxStep', 0.1, 'Events', @(t, y) stage_Separation(t, y, stages.stg1));
options_stg2 = odeset('RelTol', 1e-8, 'MaxStep', 0.1, 'Events', @(t, y) orbit_insertion(t, y));

%% First stage simulation
tspan_stg1 = 0:1/stages.stg1.u_freq:1e4;

T1 = zeros(length(tspan_stg1), 1);
Y1 = zeros(length(tspan_stg1), 6);
delta_vec_stg1 = zeros(length(tspan_stg1), 1);

Y1(1, :) = y0_stg1;
idx = 2;

alpha = 0;
theta = 0;
thetaDot = 0;


for ii = 2:length(tspan_stg1)-1

    if Y1(idx-1, 2) < 50e3
        angle = 90;
    else
        angle = 45;
    end
    err = theta - deg2rad(angle);

    delta = -stages.stg1.k1*err - stages.stg1.k2*thetaDot - stages.stg1.k3*alpha;
    if abs(delta) > stages.stg1.deltaMax 
        delta = stages.stg1.deltaMax*sign(delta);
    end
    
    [tvect, yvect, tevent, yevent,~] = ode113(@(t,y) dyn(t, y, stages.stg1, params, 1, delta), [tspan_stg1(ii) tspan_stg1(ii+1)], Y1(idx-1,:), options_stg1);
    parout = recallOdeFcn(tvect, yvect, stages.stg1, params, 1, delta);
    T1(idx:idx+length(tvect)-1) = tvect;
    Y1(idx:idx+length(tvect)-1, :) = yvect;
    delta_vec_stg1(idx:idx+length(tvect)-1) = repmat(delta, [length(tvect) 1]);
    idx = idx+length(tvect);

    alpha = parout.alpha(end);
    theta = yvect(end, 5);
    thetaDot = yvect(end, 6);

    % Exit the loop if we have stage separation
    if ~isempty(tevent)
        break
    end
end
T1(idx:end, :) = [];
Y1(idx:end, :) = [];
delta_vec_stg1(idx:end, :) = [];
clear dyn

%% Second stage simulation
tspan_stg2 = 0:1/stages.stg2.u_freq:1e4;

T2 = zeros(length(tspan_stg2), 1);
Y2 = zeros(length(tspan_stg2), 6);
delta_vec_stg2 = zeros(length(tspan_stg2), 1);

Y2(1, :) = Y1(end,:);
idx = 2;

alpha = 0;
theta = 0;
thetaDot = 0;


for ii = 2:length(tspan_stg2)-1

    if Y2(idx-1, 2) < 300e3
        angle = 45;
    else
        angle = 0;
    end
    err = theta - deg2rad(angle);

    delta = -stages.stg2.k1*err - stages.stg2.k2*thetaDot - stages.stg2.k3*alpha;
    if abs(delta) > stages.stg2.deltaMax 
        delta = stages.stg2.deltaMax*sign(delta);
    end
    
    [tvect, yvect, tevent, yevent,~] = ode113(@(t,y) dyn(t, y, stages.stg2, params, 2, delta), [tspan_stg2(ii) tspan_stg2(ii+1)], Y2(idx-1,:), options_stg2);
    parout = recallOdeFcn(tvect, yvect, stages.stg2, params, 2, delta);
    T2(idx:idx+length(tvect)-1) = tvect;
    Y2(idx:idx+length(tvect)-1, :) = yvect;
    delta_vec_stg2(idx:idx+length(tvect)-1) = repmat(delta, [length(tvect) 1]);
    idx = idx+length(tvect);

    alpha = parout.alpha(end);
    theta = yvect(end, 5);
    thetaDot = yvect(end, 6);

    % Exit the loop if we have stage separation
    if ~isempty(tevent)
        break
    end
end
T2(idx:end, :) = [];
Y2(idx:end, :) = [];
delta_vec_stg2(idx:end, :) = [];
clear dyn

%% Retrieve data from ode
clc
T = [T1; T2+T1(end)];
Y = [Y1; Y2];
delta_vec = [delta_vec_stg1; delta_vec_stg2];

qdyn = zeros(length(T), 1);
acc = zeros(length(T), 2);
alpha = zeros(length(T), 1);
moment = zeros(length(T), 1);
for ii = 1:length(T)
    if ii <= length(T1)
        [~, parout] = dyn(T(ii), Y(ii, :), stages.stg1, params, 1, delta_vec(ii));
    else
        [~, parout] = dyn(T(ii), Y(ii, :), stages.stg2, params, 2, delta_vec(ii));
    end
    qdyn(ii) = parout.qdyn;
    acc(ii,:) = parout.acc;
    % if isfield(parout, "t_turn") && ~isnan(parout.t_turn)
    %     t_turn = parout.t_turn;
    % end
    alpha(ii) = parout.alpha;
    moment(ii) = parout.moment;
%     gravity(ii) = parout.g;
     mass(ii) = parout.m;
%      %Fxz(ii,ii) = parout.Fxz;
%         F_x(1,ii)=parout.F_x; 
%         F_z(1,ii) = parout.F_z;
%         T_x(ii) = parout.T_x;
%         T_z(ii) = parout.T_z;
%         D(ii) = parout.D;
%         L(ii) = parout.L;
%          W(ii) = parout.W;
% 
%          T_M(ii) = parout.T;
%          D_M(ii) = parout.D_M;
%          L_M(ii)=parout.L_M;
%          Aer_arm(ii)=parout.Aer_arm;
% T_arm(ii)= parout.T_arm;
% A_M(ii)=parout.A_M;
% 
% 
%           Fxb(1,ii)=parout.F_x_body; 
%         Fzb(1,ii) = parout.F_z_body;
%         Txb(ii) = parout.Tx_body;
%         Tzb(ii) = parout.Tz_body;
%         Db(ii) = parout.D_body;
%         Lb(ii) = parout.L_body;
%          Wxb(ii) = parout.Wx_body;
%          Wzb(ii) = parout.Wz_body;

end 
downrange = params.Re./(params.Re+Y(:,1)) .* Y(:,1);
%% ELABORATION:
mass=mass';
% gravity=gravity';
% F_x = F_x';
% F_z = F_z';
% D=D';
% L=L';
% W=W';
% T_x=T_x';
% T_z=T_z';
% 
% Fxb = Fxb';
% Fzb = Fzb';
% Db=Db';
% Lb=Lb';
% Wxb=Wxb';
% Wzb=Wzb';
% Txb=Txb';
% Tzb=Tzb';
% 
% T_M=T_M';
% D_M=D_M';
% L_M=L_M';
% Aer_arm = Aer_arm';
% T_arm = T_arm';
% A_M=A_M';
% 
%% Structural anal
% 
mat_id_1 = 6;
mat_id_2 = 4;
[MAT.MAT_1] = material_selection(mat_id_1); % 1st stage and interstage
[MAT.MAT_2] = material_selection(mat_id_2); % 2nd stage and fairing

MASS=load('MASS.mat');
 
OUT.FS = 1.5;

%%  MAX Q:
Max_Q = max(qdyn);
FORCES.q = Max_Q;
 t_max_Q = find(Max_Q==qdyn);
 t_max_Q = t_max_Q(1,1);
 OUT.m = mass(t_max_Q);
alpha_Q = alpha(t_max_Q);

Max_M = max(moment);
t_max_M = find(Max_M==moment);
t_max_M = t_max_M(1,1);

OUT.delta=max(delta_vec);


% GEOMETRY.Diam_1 = 1.6;
% Diam_1 =GEOMETRY.Diam_1;
% GEOMETRY.Diam_2 = 1.3;
% Diam_2 =GEOMETRY.Diam_2;
alpha_Q = alpha(t_max_Q);

mp1 = stages.stg1.m_prop - abs((OUT.m - stages.stg1.m0)) - 109.9429;

[~,~,GEOMETRY]=COM_MoI_Stk1_dummy(mp1,MASS.MASS);

% GEOMETRY.b1=2; % fairing length
% GEOMETRY.b2=3; % forward skirt length
% GEOMETRY.b3=1; % height of fuel tank2
% GEOMETRY.b4=2; %  height of oxygen tank2
% GEOMETRY.b5=3; %  aft skirt length 2
% GEOMETRY.b6=3; % interstage length
% GEOMETRY.b7=2; % height of fuel tank1
% GEOMETRY.b8=2; % height of oxygen tank1
% GEOMETRY.b9=3; % aft skirt length 1
% GEOMETRY.b10=4; % h engine

% GEOMETRY.m1=50; 
% GEOMETRY.m2=50;
%  GEOMETRY.m3=27;
% GEOMETRY.m4=49;
% GEOMETRY.m5=145; 
%  GEOMETRY.m6=32;
% GEOMETRY.m7=11222;
% GEOMETRY.m8=321;
% GEOMETRY.m9=12;
% GEOMETRY.m10=12; 

% GEOMETRY.x_com1=1;
% GEOMETRY.x_com2=3;
% GEOMETRY.x_com3=5.5;
% GEOMETRY.x_com4=7.5;
% GEOMETRY.x_com5=10;
% GEOMETRY.x_com6=11;
% GEOMETRY.x_com7=12;
% GEOMETRY.x_com8=13;
% GEOMETRY.x_com9=14.5;
% GEOMETRY.x_com10=16.5;

FORCES.FS = 1.25;

%% Import valori Pietro:  h_it, M_it, M1, M2

h_it=load("h_it.mat");
M_it=load("M_it.mat");
M1=load("M1.mat");
M2=load("M2.mat");

n_choice = 34;
h_it_case = h_it.h_it(1,n_choice);
M_it_case=M_it.M_it(n_choice);

GEOMETRY = GEO(h_it_case,M_it_case,M1.M1,M2.M2);
%%

OUT.Cl_nose = stages.stg1.Cl;
OUT.Cl_fin = stages.stg1.Cl;
OUT.Cd = stages.stg1.Cd;

FORCES.T =stages.stg1.Thrust;

%q = qdyn(t_max_Q);
q=qdyn(t_max_M);
%alpha_n=0;
alpha_n=deg2rad(20);
AER=Aer_force(q,alpha_n);
AER.q =q;
OUT.FLAG = 2; % max q: 1, Pitch up is 2
%[THICKNESS]= Final_syst_case2_5(GEOMETRY,OUT,FORCES,MAT,AER);
%[THICKNESS]= Final_syst_case2(GEOMETRY,OUT,FORCES,MAT,AER);
%[THICKNESS]= Final_syst_case3(GEOMETRY,OUT,FORCES,MAT,AER); % 4 points

OUT.alpha=alpha_n;
%FORCES.nx=7;
FORCES.nx=1.1;
FORCES.nz=1.1;
%[THICKNESS]= Final_syst_case2nxz(GEOMETRY,OUT,FORCES,MAT,AER);
%[THICKNESS]= Final_syst_case2nxzT(GEOMETRY,OUT,FORCES,MAT,AER);
[THICKNESS]= Final_syst_case1nxz(GEOMETRY,OUT,FORCES,MAT,AER); % 1 point
%[THICKNESS]= Final_syst_case1nxzT(GEOMETRY,OUT,FORCES,MAT,AER); % 1 point

%% HANDLING:
g0 = 9.80665; %m/s^2
m2h =GEOMETRY.m1+GEOMETRY.m2+GEOMETRY.m3+GEOMETRY.m4+GEOMETRY.m5;
m1h =GEOMETRY.m6+GEOMETRY.m7+GEOMETRY.m8;
m3h = GEOMETRY.m9+GEOMETRY.m10;
Diam_1 = GEOMETRY.Diam_1;
Diam_2 = GEOMETRY.Diam_2;
t=M_it_case.th1.rp1;
Rad_2 = Diam_1/2;
Rad_1 = Rad_2 - t;
Cd_c = 0.1;
Cl_c = 7.8;
q = 0.5*1.1*(0.6*343)^2;
D = q*(pi*Rad_2^2)*Cd_c;
L=q*(pi*Rad_2^2)*Cl_c;
L_cone = L*2/3;
L_fin = L*1/3;
nx_c = 1;
nz_c = 4;
F1_z = m1h*g0*(nz_c+1) -L_fin;
F2_z = m2h*g0*(nz_c+1) -L_cone;
F1_x = m1h*g0*(nx_c);
F2_x = m2h*g0*(nx_c);
X_COM2 = GEOMETRY.l_tot-h_it_case.stg2.CG.tot;
X_COM1 = GEOMETRY.l_tot-h_it_case.stg1.CG.tot;
D_x_com=X_COM1-X_COM2;
l_tot =10;% Ala boeing, 48.2 ft-> 14 m ma prendo meno

syms R1 F2 R2 F1 F3 a b c d e

% Define the equations
eq1 = -F1 + R1 - F2 + R2 == 0;
eq2 = F1*a - F2*b + R2*(b + c) == 0;
eq3 = F2*c - R1*(b + c) + F1*(a + b + c) == 0;
eq4 = D_x_com-a-b == 0;
eq5 = l_tot-c-b == 0;

% Solve the system of equations
sol = solve([eq1, eq2, eq3,eq4,eq5,c>2], [R1, R2, a,b,c]);
F1 = F1_z;
F2 = F2_z;

sol=subs(sol);

P_a = D + F1_x+F2_x;
P_1=P_a;
R1 = double(subs(sol.R1));
R2 = double(subs(sol.R2));


%% PLOT:
L_tot = GEOMETRY.l_tot;
a = sol.a;
b = sol.b;
c=sol.c;
d = h_it_case.stg1.CG.tot-3;
x_com = a+b;
l=a+b+c;
z = L_tot-l;
% Segment lengths
x_segments = [X_COM2,a,b,c,d];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting

% Axial Forces (P) at segment ends
P_forces = [-D-F2_x,F1_x,0,0,0];
LOAD.P_forces = [P_a,0,0];

figure;
hold on;
% for i = 1:length(x_segments)
%     plot([x(i+1), x(i+1)], [P_a, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
% end
stairs(x, [-D, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Axial Force [N]');
title('Axial Force Diagram for Worst Case Carrier Maneuver');
xlim([0,L_tot]);
grid on;

l_cone = GEOMETRY.b1;


% Segment lengths
x_segments = [l_cone,X_COM2-l_cone,a,b,c,d];
x = [0, cumsum(x_segments)]; % Cumulative lengths for plotting
x = double(subs(x));
% Axial Forces (P) at segment ends
P_forces = [L_cone,-F2_z,-F2_z-R1,F1_z-F2_z-R1,-L_fin,-L_fin];
LOAD.P_forces = [P_a,0,0];
P_forces = double(subs(P_forces));
figure;
hold on;
% for i = 1:length(x_segments)
%     plot([x(i+1), x(i+1)], [P_a, P_forces(i)], '--', 'LineWidth', 1,'Color', 'b'); % Dashed vertical lines
% end
stairs(double(x), [0, P_forces] ,'LineWidth', 2, 'Color', 'b'); % Step-style plot
xlabel('x [m]');
ylabel('Shear Force [N]');
title('Shear Force Diagram for Worst Case Carrier Maneuver');
grid on;
xlim([0,L_tot]);





% %  HANDLING STATIC
% b = linspace(5,10,7);
% %b_vect = [GEOMETRY.b1;GEOMETRY.b2;GEOMETRY.b3;GEOMETRY.b4;GEOMETRY.b5;GEOMETRY.b6;GEOMETRY.b7;GEOMETRY.b8;GEOMETRY.b9;GEOMETRY.b10];
% MASS.m_tot = 16*10^3;
% GEOMETRY.L_tot=10; % Ala boeing, 48.2 ft-> 14 m ma prendo meno
% MAT = material_selection(3);
% 
% Diam_1 = GEOMETRY.Diam_1;
% Diam_2 = GEOMETRY.Diam_2;
% MASS.E = MAT.E;
% t=M_it_case.th1.rp1;
% Rad_2 = Diam_1/2;
% Rad_1 = Rad_2 - t;
% GEOMETRY.RLV=Rad_2;
% MASS.J = (pi/2)*(Rad_2^4 - Rad_1^4);
% 
% for i=1:length(b)
% 
% CLAMP(i)=Handling(MASS,GEOMETRY,b(i));
% 
% 
% end
% g0 = 9.80665; %m/s^2
% F=-MASS.m_tot*g0;
% l=GEOMETRY.L_tot;
% E=MASS.E;
% J=MASS.J;
% n_choice=1;
% delta_a = @(x)  ( (-F*b(n_choice).*x.^3) + F*b(n_choice)*(l^2 - b(n_choice)^2).*x )/(6*l*E*J);
% delta_b = @(x)  (F*(x-l).^3)/(6*l*E*J);
% x_vec_a = linspace(0,CLAMP(n_choice).a,1000);
% x_vec_b = linspace(b(n_choice),l,1000);
% phi_b=CLAMP(n_choice).phi_b;
% 
% 
% %%
% 
% function p = cubic_interp_with_derivatives(x1, y1, slope1, x2, y2, slope2)
%     A = [x1^3 x1^2 x1 1; 3*x1^2 2*x1 1 0; x2^3 x2^2 x2 1; 3*x2^2 2*x2 1 0];
%     b = [y1; slope1; y2; slope2];
%     coeffs = A\b;
% 
%     p = @(x) polyval(coeffs, x);
% end
% 
% x1 = x_vec_a(end);
% y1 = delta_a(x1);
% x2 = l;
% y2 = 0;
% slope2 = phi_b; % Slope at x2
% x_m=CLAMP(n_choice).a;
% slope1 =(1/(6*l*E*J)) .* ( -3*F*b(n_choice)*x_m^2 + F*b(n_choice)*(l^2 - b(n_choice)^2) );  % Slope at x2
% 
% p = cubic_interp_with_derivatives(x1, y1, slope1, x2, y2, slope2);
% 
% x_eval = linspace(x1, x2, 100);
% y_eval = p(x_eval);
% 
% y_min=min(y_eval);
% position_min = find(y_eval==y_min);
% x_min = x_eval(position_min);
% 
% loads.nx=2;
% loads.nz=2;
% TR.a=CLAMP(n_choice).a;
% TR.b=b(n_choice);
% TR.l=l;
% TR.J=J;
% TR.Diam1=Diam_1;
% TR.m_tot=MASS.m_tot;
% TR.x_min=x_min;
% TR.t = M_it_case.th1.rp1;
% TR.FS=OUT.FS;
% [LOAD,TR] = Attachments(TR,loads);
% 
% figure()
% plot(x_vec_a,delta_a(x_vec_a),'Color','b');
% hold on;
% plot(x_eval, y_eval,'Color','b');
% hold on;
% plot([0, l], [-l/3000, -l/3000], 'Color', 'r', 'LineStyle', '--');
% plot(x_min,y_min,'o','MarkerEdgeColor','k');
% xlabel('x [m]');
% ylabel('Elastic deflection [m]');
% title('Elastic deflection of the Launch Vehicle in hanging configuration');
% legend('','Elastic Deflection','Deflection limit','Baricenter');
% grid on;
% 
% SOL=CLAMP(n_choice);
% 
% x_com =GEOMETRY.x_com;
% SOL.L_clamp = SOL.A/(2*pi*Rad_2);
% 
% SOL.x_max_LV=x_com;
% 
% SOL.x_a = x_com -SOL.a;
% SOL.x_b = x_com +SOL.b;


% 
% %% MAX Q:
% Max_Q = max(qdyn);
% q = Max_Q;
% t_max_Q = find(Max_Q==qdyn);
% t_max_Q = t_max_Q(1,1);
% 
% % OUT1.rot_angle= rot_angle(t_max_Q);
% % OUT1.D= D(t_max_Q);
% % OUT1.gamma= gamma(t_max_Q);
% % OUT1.L=L(t_max_Q);
% % OUT1.delta=delta_vec(t_max_Q);
% % OUT1.theta=theta(t_max_Q);
% % OUT1.g =gravity(t_max_Q);
% % OUT1.m=mass(t_max_Q);
% % OUT1.T = Thrust(t_max_Q);
% 
% OUT1.T_M=T_M(t_max_Q)*sin(delta_vec(t_max_Q));
% OUT1.D_M=D_M(t_max_Q)*sin(alpha(t_max_Q));
% OUT1.L_M=L_M(t_max_Q)*cos(alpha(t_max_Q));
% OUT1.Aer_arm=Aer_arm(t_max_Q);
% OUT1.T_arm=T_arm(t_max_Q);
% OUT1.A_M=A_M(t_max_Q);
% 
% OUT1.D = D(t_max_Q);
% OUT1.L=L(t_max_Q);
% OUT1.T_x= T_x(t_max_Q);
% OUT1.T_z= T_z(t_max_Q);
% OUT1.F_x= F_x(t_max_Q);
% OUT1.F_z= F_z(t_max_Q);
% OUT1.W= W(t_max_Q);
% OUT1.m = mass(t_max_Q);
% 
% OUT1.Db = Db(t_max_Q);
% OUT1.Lb=Lb(t_max_Q);
% OUT1.Txb= Txb(t_max_Q);
% OUT1.Tzb= Tzb(t_max_Q);
% OUT1.Fxb= Fxb(t_max_Q);
% OUT1.Fzb= Fzb(t_max_Q);
% OUT1.Wxb= Wxb(t_max_Q);
% OUT1.Wzb= Wzb(t_max_Q);
% OUT1.m = mass(t_max_Q);
% 
% OUT1.D = Db(t_max_Q);
% OUT1.L=Lb(t_max_Q);
% OUT1.T_x= Txb(t_max_Q);
% OUT1.T_z= Tzb(t_max_Q);
% OUT1.F_x= Fxb(t_max_Q);
% OUT1.F_z= Fzb(t_max_Q);
% OUT1.W_x= Wxb(t_max_Q);
% OUT1.W_z= Wzb(t_max_Q);
% 
% OUT1.moment=moment(t_max_Q);
% 
% M_balance = OUT1.moment - (OUT1.T_M*OUT1.T_arm) - (OUT1.D_M*OUT1.Aer_arm)+ (OUT1.L_M*OUT1.Aer_arm);
% M_balance2 = OUT1.moment - (OUT1.T_M*OUT1.T_arm) - (OUT1.A_M*OUT1.Aer_arm);
% 
% % X_balance = OUT1.F_x - OUT1.T_x -OUT1.D;
% % 
% % X_b_balance = OUT1.Fxb - OUT1.Txb - OUT1.Db - OUT1.Wxb;
% % 
% % Z_balance = OUT1.F_z - OUT1.T_z -OUT1.L -OUT1.W;
% % 
% % Z_b_balance = OUT1.Fzb - OUT1.Tzb - OUT1.Lb - OUT1.Wzb;
% 
% %%
% mp1 = stages.stg1.m_prop - abs((OUT1.m - stages.stg1.m0)) - 109.9429;
% 
% [~,~,GEOMETRY1]=COM_MoI_Stk1_dummy(mp1,MASS.MASS);
% 
% % GEOMETRY1.b1 = 2;
% % GEOMETRY1.b2 = 3;
% % GEOMETRY1.b3 = 3;
% % GEOMETRY1.b4 = 1;
% % GEOMETRY1.x_com =8;
% 
% %[BALANCE]= Load_Distr(OUT1,GEOMETRY1);
% [BALANCE]= Load_Distr_BODY(OUT1,GEOMETRY1);
% 
% %%
% 
% FORCES1.F_x = F_x(t_max_Q);
% FORCES1.F_z = F_z(t_max_Q);
% 
% P_balance1 = FORCES1.F_x- FORCES1.D_x - FORCES1.L_x - FORCES1.T_x -FORCES1.G_x;
% R_balance1 = FORCES1.F_z- FORCES1.D_z - FORCES1.L_z - FORCES1.T_z -FORCES1.G_z;
% M_1 = FORCES1.M;
% 
% [THICK1]= Stran_Pitchup(FORCES1,GEOMETRY,MAT,OUT1)
% 
% %check_1 =GEOMETRY.m1+GEOMETRY.m2+GEOMETRY.m3+GEOMETRY.m4 - OUT1.m
% 
% %% Pitch-up Maneuver: (mp2>mp1)!!!!
% Max_M = max(moment);
% t_max_M = find(Max_M==moment);
% t_max_M = t_max_M(1,1);
% 
% 
% OUT2.rot_angle= rot_angle(t_max_M);
% OUT2.D= D(t_max_M);
% OUT2.gamma= gamma(t_max_M);
% OUT2.L=L(t_max_M);
% OUT2.delta=delta_vec(t_max_M);
% OUT2.theta=theta(t_max_M);
% OUT2.g=gravity(t_max_M);
% OUT2.m=mass(t_max_M);
% OUT2.T = Thrust(t_max_M);
% OUT2.moment=moment(t_max_M);
% 
% mp2 = stages.stg1.m_prop - abs((OUT2.m - stages.stg1.m0)) - 109.9429;
% 
% [~,~,GEOMETRY]=COM_MoI_Stk1_dummy(mp2,MASS.MASS);
% 
% 
% FORCES2 = Body_frame_Forces(OUT2);
% 
% F_x2 = F_x(t_max_M);
% F_z2 = F_z(t_max_M);
% 
% P_balance2 = F_x2- FORCES2.D_x - FORCES2.L_x - FORCES2.T_x -FORCES2.G_x;
% R_balance2 = F_z2- FORCES2.D_z - FORCES2.L_z - FORCES2.T_z -FORCES2.G_z;
% M_2 = FORCES2.M;
% 
% chech2=GEOMETRY.m1+GEOMETRY.m2+GEOMETRY.m3+GEOMETRY.m4 - OUT2.m




%% Event functions

function [value, isterminal, direction] = stage_Separation(t, ~, stage)
    value = t - (stage.t_burn_tot + stage.t_wait + 1);
    isterminal = 1;
    direction = 1;
end

function [value, isterminal, direction] = orbit_insertion(~, y)
    value = y(4);
    isterminal = 1;
    direction = 0;
end

