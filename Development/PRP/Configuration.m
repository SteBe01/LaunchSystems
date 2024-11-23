clear
close all
clc

% INPUT
d_rocket=1.25; %m
d_old=0.25;
N=7;
%L=0.67; old
L_engine=0.8;
gimbal = 5; %degs

%COMPUTATION
safety = L_engine*sind(gimbal);
angle_poly = 90*(1-2/N);
d_new = (d_rocket*cosd(angle_poly)-safety)/(1+cosd(angle_poly));

% PLOT
figure
xx = 0:pi/180:2*pi;
plot(d_rocket/2*cos(xx), d_rocket/2*sin(xx),'b', LineWidth=3)
hold on
axis equal
pos_center_x = [ 0 (d_rocket-d_new)/2*cos(linspace(0,2*pi,N+1))];
pos_center_y = [ 0 (d_rocket-d_new)/2*sin(linspace(0,2*pi,N+1))];
for i=1:N+1
    plot( d_old/2        *cos(xx) + pos_center_x(i),  d_old/2        *sin(xx) + pos_center_y(i), 'g')
    plot( d_new/2        *cos(xx) + pos_center_x(i),  d_new/2        *sin(xx) + pos_center_y(i), 'k', LineWidth=3)
    plot((d_new/2+safety)*cos(xx) + pos_center_x(i), (d_new/2+safety)*sin(xx) + pos_center_y(i), '--r')
end
legend("Rocket contour", "Old engine", "New engine", "Safety radius")
title(strcat("Engine configuration: ", num2str(N), " engines in the outer ring, total of ", num2str(N+1)))
eps_old=14;
eps_new=d_new^2/d_old^2*eps_old;


