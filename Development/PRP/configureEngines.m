function eps_new = configureEngines (d_rocket, N_engines, L_engine, gimbal, d_th, flag)

% INPUT
N_out=N_engines-1;

%COMPUTATION
safety = L_engine*sind(gimbal);
angle_poly = 90*(1-2/N_out);
d_new = (d_rocket*cosd(angle_poly)-safety)/(1+cosd(angle_poly));

% PLOT
if flag
    figure
    xx = 0:pi/180:2*pi;
    plot(d_rocket/2*cos(xx), d_rocket/2*sin(xx),'b', LineWidth=3)
    hold on
    axis equal
    pos_center_x = [ 0 (d_rocket-d_new)/2*cos(linspace(0,2*pi,N_out+1))];
    pos_center_y = [ 0 (d_rocket-d_new)/2*sin(linspace(0,2*pi,N_out+1))];
    for i=1:N_out+1
        plot( d_new/2        *cos(xx) + pos_center_x(i),  d_new/2        *sin(xx) + pos_center_y(i), 'k', LineWidth=3)
        plot((d_new/2+safety)*cos(xx) + pos_center_x(i), (d_new/2+safety)*sin(xx) + pos_center_y(i), '--r')
    end
    legend("Rocket contour", "New engine", "Safety radius")
    title(strcat("Engine configuration: ", num2str(N_out), " engines in the outer ring, total of ", num2str(N_out+1)))
end

eps_new=d_new^2/d_th^2;


