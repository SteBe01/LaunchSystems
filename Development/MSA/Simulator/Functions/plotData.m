function plotData(T, Y, params, parout, idxStage)
    set(0, 'DefaultLineLineWidth', 1.5)

    T1 = T(1:idxStage);
    Y1 = Y(1:idxStage, :);

    g_vec = 398600*1e9./vecnorm(Y(:,1:2),2,2).^2;
    downrange = params.Re./(params.Re+Y(:,1)) .* Y(:,1);

    beta = atan2(Y(:,2), Y(:,1));
    ang = pi/2-beta;
    vec_rotated = zeros(length(T),2);
    for ii = 1:length(T)
        dcm = [cos(ang(ii)) -sin(ang(ii)); sin(ang(ii)) cos(ang(ii))];
        vec_rotated(ii, :) = dcm*Y(ii, 3:4)';
    end

    boundary = 0;

    subplot(2,2,1), hold on, grid on, title("Radius over time"), xlabel("Time [s]"), ylabel("Radius [km]")
    plot(T, (vecnorm(Y(:,1:2),2,2) - params.Re)/1e3)
    xline(T1(end), '--k', 'Staging')
    subplot(2,2,2), hold on, grid on, title("Downrange over time"), xlabel("Time [s]"), ylabel("Downrange [km]")
    plot(T, downrange/1e3)
    xline(T1(end), '--k', 'Staging')
    subplot(2,2,3), hold on, grid on, title("Horizontal (rotated) velocity over time"), xlabel("Time [s]"), ylabel("Velocity [km/s]")
    plot(T, vec_rotated(1:end-boundary,1)/1e3)
    xline(T1(end), '--k', 'Staging')
    subplot(2,2,4), hold on, grid on, title("Vertical (rotated) velocity over time"), xlabel("Time [s]"), ylabel("Velocity [km/s]")
    plot(T, vec_rotated(1:end-boundary,2)/1e3)
    xline(T1(end), '--k', 'Staging')

    figure, hold on, grid on, title("Dynamic pressure wrt altitude for first stage"), xlabel("Altitude [km]"), ylabel("Qdyn [kPa]")
    plot(Y1(1:end, 2)/1e3 - params.Re/1e3, parout.qdyn(1:parout.idxStg1)/1e3);

    figure
    subplot(4,1,1), hold on, grid on, title("$\xi$ evolution", 'Interpreter', 'latex'), xlabel("Time [s]"), ylabel("Xi [deg]")
    plot(T, unwrap(rad2deg(Y(:,5)) + 90 - rad2deg(atan2(Y(:,2),Y(:,1)))))
    xline(T1(end), '--k', 'Staging')
    subplot(4,1,2), hold on, grid on, title("$\dot\theta$ over time", 'Interpreter', 'latex'), xlabel("Time [s]"), ylabel("Theta dot [deg/s]")
    plot(T, rad2deg(Y(:, 6)))
    xline(T1(end), '--k', 'Staging')
    subplot(4,1,3), hold on, grid on, title("$\alpha$ evolution", 'Interpreter', 'latex'), xlabel("Time [s]"), ylabel("Alpha [deg]")
    plot(T, unwrap(rad2deg(parout.alpha)))
    xline(T1(end), '--k', 'Staging')
    subplot(4,1,4), hold on, grid on, title("$\delta$ evolution", 'Interpreter', 'latex'), xlabel("Time [s]"), ylabel("Delta [deg]")
    plot(T, rad2deg(parout.delta_vec))
    xline(T1(end), '--k', 'Staging')

    figure
    subplot(2,2,1), hold on, grid on, title("Axial acceleration over time"), xlabel("Time [s]"), ylabel("Acceleration [g]")
    plot(T, parout.acc(:,1)/params.g0)
    xline(T1(end), '--k', 'Staging')
    subplot(2,2,2), hold on, grid on, title("Normal acceleration over time"), xlabel("Time [s]"), ylabel("Acceleration [g]")
    plot(T, parout.acc(:,2)/params.g0)
    xline(T1(end), '--k', 'Staging')
    subplot(2,2,3), hold on, grid on, title("Acceleration norm over time"), xlabel("Time [s]"), ylabel("Acceleration norm [g]")
    plot(T, vecnorm(parout.acc,2,2)./params.g0)
    plot(T, g_vec./params.g0, '--')
    xline(T1(end), '--k', 'Staging')
    subplot(2,2,4), hold on, grid on, title("Moments over time"), xlabel("Time [s]"), ylabel("Moment [Nm]")
    plot(T, parout.moment)
    xline(T1(end), '--k', 'Staging')

    figure, hold on, grid on, title("Propellant mass over time"), xlabel("Time [s]"), ylabel("Fuel mass [kg]")
    plot(T, Y(:,end))
    xline(T1(end), '--k', 'Staging')

end

