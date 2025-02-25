%% Animation

factor = 5;

clc

len = length(parout.dcm);

ii_new = 1:factor:len;

F_L_in = parout.F_L_in;
norma_old = 1;
max_L = 1;
for ii = ii_new
    F_L_in_temp = F_L_in(ii,:);
    norma = norm(F_L_in_temp);
    if norma > norma_old
        max_L = norma;
        norma_old = max_L;
    end
end
F_D_in = parout.F_D_in;
norma_old = 1;
max_D = 1;
for ii = ii_new
    F_D_in_temp = F_D_in(ii,:);
    norma = norm(F_D_in_temp);
    if norma > norma_old
        max_D = norma;
        norma_old = max_D;
    end
end


figure, hold on, grid on, axis equal, xlim([-1 1]), ylim([-1 1])
for ii = ii_new
    dcm = parout.dcm;
    dcm = dcm(:,:,ii);
    F_L_in = parout.F_L_in;
    F_D_in = parout.F_D_in;
    F_L_in = F_L_in(ii,:);
    F_D_in = F_D_in(ii,:);
    F_L_in = F_L_in./max_L;
    F_D_in = F_D_in./max_D;

    plt1 = plot([dcm(1,1) 0],[dcm(2,1) 0], Color="black");  % body axis
    plt2 = plot([dcm(1,2) 0],[dcm(2,2) 0], Color="black");  % body axis

    wind_norm = norm(Y(ii,3:4));
    plt3 = plot([Y(ii,3)./wind_norm 0],[Y(ii,4)./wind_norm 0], Color="blue");  % wind

    plt4 = plot([F_L_in(1) 0],[F_L_in(2) 0], Color="red");
    plt5 = plot([F_D_in(1) 0],[F_D_in(2) 0], Color="red");

    drawnow
    fprintf("Theta: %3.1f deg\n", rad2deg(Y(ii, 5)))
    fprintf("%d, Position: vz = %1.3f km/s, vx = %1.3f km/s, h = %3.1f km\n", T(ii), Y(ii, 4)./1e3, Y(ii, 3)./1e3, Y(ii, 2)./1e3)

    % pause(0.01)
    
    delete(plt1)
    delete(plt2)
    delete(plt3)
    delete(plt4)
    delete(plt5)
end

close all

