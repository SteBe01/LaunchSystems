
%% Plot orbits
figure, hold on, grid on, title("Radius [km]"), xlabel("Time [s]"), ylabel("Radius [km]")
for ii = 1:N_sim
    plot(saveSim{ii}.t, (vecnorm(saveSim{ii}.Y(:,1:2),2,2)-nom_params.Re)/1e3)
end

%% Plot Downrange

figure, hold on, grid on, title("Downrange [km]"), xlabel("Time [s]"), ylabel("Downrange [km]")
for ii = 1:N_sim
    downrange = nom_params.Re./(nom_params.Re+saveSim{ii}.Y(:,1)) .* saveSim{ii}.Y(:,1);
    plot(saveSim{ii}.t, downrange/1e3)
end


%% Orbit radius error
max_alt = zeros(N_sim,1);
mu_alt = zeros(N_sim,1);
std_alt = zeros(N_sim,1);

for ii = 1:N_sim
    max_alt(ii) = norm(saveSim{ii}.Y(end,1:2));
    mu_alt(ii) = mean(max_alt(1:ii));
    std_alt(ii) = std(max_alt(1:ii));
end
figure; hold on; grid on;
subplot(2,1,1);
plot(1:N_sim, (mu_alt-6778e3)/1e3);
yline((mean(max_alt)-6778e3)/1e3);
title("Mean orbit radius error"); ylabel("Error [km]"); xlabel("Number of simulations")
subplot(2,1,2);
plot(1:N_sim, std_alt/1e3);
yline(std(max_alt)/1e3);
title("std orbit radius error"); ylabel("Error std [km]"); xlabel("Number of simulations")

%% Velocity error

max_vel = zeros(N_sim,1);
mu_vel = zeros(N_sim,1);
std_vel = zeros(N_sim,1);

for ii = 1:N_sim
    max_vel(ii) = norm(saveSim{ii}.Y(end,3:4));
    mu_vel(ii) = mean(max_vel(1:ii));
    std_vel(ii) = std(max_vel(1:ii));
end
v_orb = sqrt(398600/6778);
figure; hold on; grid on;
subplot(2,1,1);
plot(1:N_sim, (mu_vel/1e3-v_orb));
yline((mean(max_vel)/1e3-v_orb));
title("Mean orbit velocity error"); ylabel("Error [km/s]"); xlabel("Number of simulations")
subplot(2,1,2);
plot(1:N_sim, std_vel/1e3);
yline(std(max_vel)/1e3);
title("std orbit velocity error"); ylabel("Error std [km/s]"); xlabel("Number of simulations")

%% Plot downrange over geoplot

inclination = deg2rad(98);
ang = inclination - pi/2;
Rot_matrix = [cos(ang) -sin(ang); sin(ang) cos(ang)];
lat0 = 34.751330;
lon0 = -120.52023;
h0 = 112.1;

figure;

downrange_end = zeros(N_sim, 1);
coords_end = zeros(N_sim, 2);

for ii = 1:N_sim
    n = length(saveSim{ii}.t);
    Rot_mat = repmat(Rot_matrix, 1,1,n);
    downrange = -(nom_params.Re./(nom_params.Re+saveSim{ii}.Y(:,1)) .* saveSim{ii}.Y(:,1) + 500e3);
    pos = pagemtimes(Rot_mat, reshape([downrange zeros(n,1)]', 2, 1, n));
    pos = squeeze(pos)';
    downrange_end(ii) = downrange(end);

    coords = zeros(n, 2);
    [coords(:,1), coords(:,2)] = ned2geodetic(pos(:,1), pos(:,2), h0*ones(n, 1), lat0, lon0, h0, wgs84Ellipsoid);

    % geoplot(coords(:,1), coords(:,2), 'LineStyle', '-');
    geoplot(coords(end,1), coords(end,2), 'LineStyle','none', 'Marker','.', 'MarkerSize',12);
    hold on;
    coords_end(ii, :) = coords(end,1:2);
end
geobasemap satellite

coord_mean = mean(coords_end);
coord_std = std(coords_end);
geoplot(coord_mean(1), coord_mean(2), 'LineStyle','none', 'Marker','x', 'MarkerSize', 24);
% semimajor = abs(max(downrange_end) - min(downrange_end))/1e5;
% ecc = axes2ecc(semimajor,semimajor/3);
% [lat1,lon1] = ellipse1(coord_mean(1),coord_mean(2),[semimajor ecc], 8);
% geoplot(lat1,lon1,"LineWidth",2)
% geobasemap satellite