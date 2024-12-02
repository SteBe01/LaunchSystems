function plotData(T, Y, params, parout, idxStage)
    set(0, 'DefaultLineLineWidth', 1.5)

    T1 = T(1:idxStage);
    Y1 = Y(1:idxStage, :);

    g_vec = params.g0./((1+Y1(:,2)/params.Re).^2);
    downrange = params.Re./(params.Re+Y(:,1)) .* Y(:,1);


    boundary = 0;
    subplot(2,2,1), hold on, grid on, title("Downrange over time"), xlabel("Time [s]"), ylabel("Downrange [km]")
    plot(T, downrange/1e3)
    xline(T1(end), '--k', 'Staging')
    subplot(2,2,2), hold on, grid on, title("Vertical position over time"), xlabel("Time [s]"), ylabel("Altitude [km]")
    plot(T, Y(1:end-boundary,2)/1e3)
    xline(T1(end), '--k', 'Staging')
    subplot(2,2,3), hold on, grid on, title("Horizontal velocity over time"), xlabel("Time [s]"), ylabel("Velocity [km/s]")
    plot(T, Y(1:end-boundary,3)/1e3)
    xline(T1(end), '--k', 'Staging')
    subplot(2,2,4), hold on, grid on, title("Vertical velocity over time"), xlabel("Time [s]"), ylabel("Velocity [km/s]")
    plot(T, Y(1:end-boundary,4)/1e3)
    xline(T1(end), '--k', 'Staging')

    figure, hold on, grid on, title("Dynamic pressure wrt altitude"), xlabel("Altitude [km]"), ylabel("Qdyn [kPa]")
    plot(Y(1:end, 2)/1e3, parout.qdyn/1e3);
    xline(Y1(end, 1)/1e3, '--k', 'Staging')

    figure, hold on, grid on, title("Altitude wrt Downrange"), xlabel("Downrange [km]"), ylabel("Altitude [km]")
    plot(downrange/1e3, Y(1:end, 2)/1e3)
    xline(downrange(length(T1))/1e3, '--k', 'Staging')

    figure
    subplot(4,1,1), hold on, grid on, title("Theta over time"), xlabel("Time [s]"), ylabel("Theta [deg]")
    plot(T, rad2deg(Y(:, 5)))
    xline(T1(end), '--k', 'Staging')
    subplot(4,1,2), hold on, grid on, title("Theta dot over time"), xlabel("Time [s]"), ylabel("Theta dot [deg/s]")
    plot(T, rad2deg(Y(:, 6)))
    xline(T1(end), '--k', 'Staging')
    subplot(4,1,3), hold on, grid on, title("$\alpha$ evolution", 'Interpreter', 'latex'), xlabel("Time [s]"), ylabel("Alpha [deg]")
    plot(T, rad2deg(parout.alpha))
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
    plot(T, 1./((1+Y(:,2)./params.Re).^2), '--')
    xline(T1(end), '--k', 'Staging')
    subplot(2,2,4), hold on, grid on, title("Moments over time"), xlabel("Time [s]"), ylabel("Moment [Nm]")
    plot(T, parout.moment)
    xline(T1(end), '--k', 'Staging')

    figure, hold on, grid on, title("Propellant mass over time"), xlabel("Time [s]"), ylabel("Fuel mass [kg]")
    plot(T, Y(:,end))
    xline(T1(end), '--k', 'Staging')
end