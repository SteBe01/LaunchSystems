
%% run this with the relative .mat file
clc;

r = 34; %row of the data selected

%fundamental parameters recovery
AR = h_it(r).stg1.R_lox / h_it(r).stg1.dome_lox; %[-] aspect ratio of the tanks
diam1 = M_it(r).diam1; %[m] first stage diameter
r1 = diam1/2; %[m] first stage radius
diam2 = M_it(r).diam2; %[m] second stage diameter
r2 = diam2/2; %[m] second stage radius

%height recovery
%first stage:
h0 = 0; %[m] starting height 
h1 = h_it(r).stg1.motor; %[m] height of the first stage motor
h2 = h_it(r).stg1.C3 + h1; %[m] height at which the cylindrical part of the lox tank starts (first stage)
h3 = h2 + h_it(r).stg1.cyl_lox +...
    h_it(r).stg1.C2 + h_it(r).stg1.cyl_rp1; %[m] height at which the cylindrical part of the rp1 tank ends (first stage)
%h4 = h_it(r).stg1.tot; %[m] height of the first stage end
h4 = h3 + h_it(r).stg1.C1; %[m] height of the first stage end
%second stage:
h5 = h4 + h_it(r).stg2.C3; %[m] height at which the cylindrical part of the lox tank starts (second stage)
h6 = h5 + h_it(r).stg2.cyl_lox +...
    h_it(r).stg2.C2 + h_it(r).stg2.cyl_rp1; %[m] height at which the cylindrical part of the rp1 tank ends (second stage)
h7 = h6 + h_it(r).stg2.C1; %[m] height of the second stage end
%fairing:
% h8 = h7 + h_it(r).fairing - 2 * diam2; %[m] height at which the fairing cone starts
h8 = h7; %[m] height at which the fairing cone starts
h9 = h7 + h_it(r).fairing; %[m] total height of the rocket
%CoM:
hCG = h_it(r).CG; %[m] total rocket CoM position (in "h" coordinates)

%silhouette parameters:
h_rocket = [h1, h1, h2, h3, h4, h5, h6, h7, h8, h9]';
y = [ 0, r1, r1, r1, r2, r2, r2, r2, r2,  0]'; %[m]
h_rocket = [h_rocket; flip(h_rocket)]; %[m] change coordinate system to adequate to the convention
x = h9 - h_rocket; %[m]
y = [y; -flip(y)]; %[m]
xCG = h9 - hCG; %[m] total rocket CoM position (in "x" coordinates)
yCG = 0; %[m]

%silhouette plot:
figure(5)
plot(x, y, '-k', LineWidth=1); grid on; axis equal; hold on;
% plot([diam1/2, -diam1/2], [h1.til_tank-h1.dome_rp1,h1.til_tank-h1.dome_rp1], '--k');
% plot([diam2/2, -diam2/2], [h1.attach,h1.attach], '--k');
% plot([diam2/2, -diam2/2], [h3.tot-2*diam2,h3.tot-2*diam2], '--k');
plot(xCG, yCG, '+r'); 

%tank plot:
cap1 = @(k) sqrt(r1^2 - k.^2)/AR;
y1 = linspace(-r1, r1, 1e4);
x1 = @(k) h9 - ( h2 - cap1(k) );
plot(x1(y1), y1, '-k');
y2 = y1;
x2 = @(k) h9 - ( h2 + h_it(r).stg1.cyl_lox + cap1(k) );
plot(x2(y2), y2, '-k');
AR_1 = ( r1 + h_it(r).stg1.C2 ) / ( h_it(r).stg1.dome_lox + h_it(r).stg1.C2);
cap1_AR_1 = @(k) sqrt( ( r1 + h_it(r).stg1.C2 )^2 - k.^2)/AR_1;
y3 = y1;%linspace(-r1 - h_it(r).stg1.C2, r1 + h_it(r).stg1.C2, 1e4); 
x3 = @(k) h9 - ( h2 + h_it(r).stg1.cyl_lox + cap1_AR_1(k) );
plot(x3(y3), y3, '-k');
y4 = y1;
x4 = @(k) h9 - ( h3 + cap1(k) );
plot(x4(y4), y4, '-k');
cap2 = @(k) sqrt(r2^2 - k.^2)/AR;
y5 = linspace(-r2, r2, 1e4);
x5 = @(k) h9 - ( h5 - cap2(k) );
plot(x5(y5), y5, '-k');
y6 = y5;
x6 = @(k) h9 - ( h5 + h_it(r).stg2.cyl_lox + cap2(k) );
plot(x6(y6), y6, '-k');
AR_2 = ( r2 + h_it(r).stg2.C2 ) / ( h_it(r).stg2.dome_lox + h_it(r).stg2.C2);
cap2_AR_2 = @(k) sqrt( ( r2 + h_it(r).stg2.C2 )^2 - k.^2)/AR_2;
y7 = y5;%linspace(-r2 - h_it(r).stg2.C2, r2 + h_it(r).stg2.C2, 1e4); 
x7 = @(k) h9 - ( h5 + h_it(r).stg2.cyl_lox + cap2_AR_2(k) );
plot(x7(y7), y7, '-k');
y8 = y5;
x8 = @(k) h9 - ( h6 + cap2(k) );
plot(x8(y8), y8, '-k');

%fairing plot:
x9 = h9 - h8 * ones(2,1);
y9 = [-r2; r2];
plot(x9, y9, '--k',LineWidth=1);
x10 = h9 - h7 * ones(2,1);
y10 = y9;
plot(x10, y10, '--k',LineWidth=1);

%motor #i:
xm = [0, 0, (h1-0.26)*0.5, h1-0.26, h1-0.19, h1, h1, h1-0.19, h1-0.26, (h1-0.26)*0.5, 0, 0];
h1_2 = h_it(r).stg2.motor;
xm_2 = [0, 0, (h1_2-0.26)*0.5, h1_2-0.26, h1_2-0.19, h1_2, h1_2, h1_2-0.19, h1_2-0.26, (h1_2-0.26)*0.5, 0, 0];
ym = [0, 0.15, 0.12, 0.035, 0.063, 0.063, -0.063, -0.063, -0.035, -0.12, -0.15, 0];
grayColor = [128 128 128]/255;

%first stage:
alpha = linspace(0, 2*pi, M_it(r).n_mot1 + 1);
for i = 1 : M_it(r).n_mot1
    xm1 = h9 - xm;
    ym1 = (r1 - 0.15)*cos(alpha(i)) - ym;
    plot(xm1, ym1, '-','Color', grayColor, LineWidth=1);
end

%second stage:
xm2 = h9 - ( h4 + xm_2 - h1_2 );
ym2 = ym;
plot(xm2, ym2, '-','Color', grayColor, LineWidth=1);
ystg2 = y9;
xstg2 = h9 - h_it(r).stg1.tot * ones(2,1);
plot(xstg2, ystg2, '-k',LineWidth=1);
 
%adapter:
h_ad = r2 / 2;
r1_ad = r2 / 1.6;
r2_ad = r2 / 1.2;
X_a = (h9 - h7) - [0, 0, h_ad, h_ad, 0, 0];
Y_a = [0, r2_ad, r1_ad, -r1_ad, -r2_ad, 0];
plot(X_a, Y_a, '--r'); %adapter

%payload
r_pl = 0.90 * r1_ad;
theta = linspace(0, 2*pi, 3e1 + 1);
ypl = r_pl * sin(theta);
xpl = (h9 - h7 - h_ad * 0.3 - r_pl ) + r_pl * cos(theta);
plot( xpl, ypl, '--b');


