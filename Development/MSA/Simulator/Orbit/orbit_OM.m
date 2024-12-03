clc; clearvars; close all

T = 2*pi*sqrt(6778^3/398600);
[T, Y] = ode113(@(t, y) odeFun(t, y), [0 T], [0 6778 sqrt(398600/(6778)) 0], odeset('MaxStep', 1));
figure;
plot(Y(:,1), Y(:,2));
axis equal

figure;
plot(T, Y(:, 3:4));
legend("Vel_x", "Vel_z")
title("Velocities")

for ii = 1:length(T)
dY(ii,:) = odeFun(T(ii), Y(ii, :))';
end
figure;
plot(T, dY(:,3:4)*1e3);
legend("Acc_x", "Acc_z")
title("Accelerations")

figure;
plot(T, unwrap(atan2(Y(:,2), Y(:,1))))
title("$\beta$", 'Interpreter','latex')


function dY = odeFun(t,Y)
    r = Y(1:2);
    v = Y(3:4);

    dY(1:2) = v;
    dY(3:4) = -398600/norm(r)^3 * r;
    dY = dY';
end