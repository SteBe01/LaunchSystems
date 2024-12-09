
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
