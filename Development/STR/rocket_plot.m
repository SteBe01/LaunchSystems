
%% run this with the relative .mat file


r = 28; %row of the data selected

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
h4 = h_it(r).stg1.tot; %[m] height of the first stage end
%second stage:
h5 = h4 + h_it(r).stg2.C3; %[m] height at which the cylindrical part of the lox tank starts (second stage)
h6 = h5 + h_it(r).stg2.cyl_lox +...
    h_it(r).stg2.C2 + h_it(r).stg2.cyl_rp1; %[m] height at which the cylindrical part of the rp1 tank ends (second stage)
h7 = h6 + h_it(r).stg2.C1; %[m] height of the second stage end
%fairing:
h8 = h7 + h_it(r).fairing - 2 * diam2; %[m] height at which the fairing cone starts
h9 = h7 + h_it(r).fairing; %[m] total height of the rocket
%CoM:
hCG = h_it(r).CG; %[m] total rocket CoM position (in "h" coordinates)

%silhouette parameters:
h_rocket = [h0, h0, h1, h2, h3, h4, h5, h6, h7, h8, h9]';
y = [ 0, r1, r1, r1, r1, r2, r2, r2, r2, r2,  0]'; %[m]
h_rocket = [h_rocket; flip(h_rocket)]; %[m] change coordinate system to adequate to the convention
x = h9 - h_rocket; %[m]
y = [y; -flip(y)]; %[m]
xCG = h9 - hCG; %[m] total rocket CoM position (in "x" coordinates)
yCG = 0; %[m]

%silhouette plot:
figure(4)
plot(x, y, '-k'); grid on; axis equal; hold on;
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
cap3 = @(k) sqrt( ( r1 + h_it(r).stg1.C2 )^2 - k.^2)/AR_1;
y3 = y1;
x3 = @(k) h9 - ( h2 + h_it(r).stg1.cyl_lox + cap3(k) );
plot(x3(y3), y3, '-k');



