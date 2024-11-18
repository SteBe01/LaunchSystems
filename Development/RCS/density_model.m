function rho = density_model(h)

h = h/1000; %[km]

% if h >= 0 && h < 25
%     h0 = 0;
%     rho0 = 1.225;
%     H = 7.249;
% 20
% elseif h >= 25 && h < 30
%     h0 = 25;
%     rho0 = 3.899e-2;
%     H = 6.349;
%    
% elseif h >= 30 && h < 40
%     h0 = 30;
%     rho0 = 1.774e-2;
%     H = 6.682;
% elseif h >= 40 && h < 50
%     h0 = 40;
%     rho0 = 3.972e-3;
%     H = 7.554;
% elseif h >= 50 && h < 60
%     h0 = 50;
%     rho0 = 1.057e-3;
%     H = 8.4;
% elseif h >= 60 && h < 70
%     h0 = 60;
%     rho0 = 3.206e-4;
%     H = 7.714;
% elseif h >= 70 && h < 80
%     h0 = 70;
%     rho0 = 8.770e-5;
%     H = 6.549;
% elseif h >= 80 && h < 90
%     h0 = 80;
%     rho0 = 1.905e-5;
%     H = 5.799;
%     elseif h >= 90 
%     h0 = 90;
%     rho0 = 3.396e-6;
%     H = 5.382;
% end
h0 = 0;
    rho0 = 1.225;
    H = 7.249;
rho = rho0* exp(-(h- h0)/H);



end