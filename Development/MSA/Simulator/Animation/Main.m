%% Animation

clc
close all

clearvars("TR_engine_vect", "TR_body_vect")

theta_vect = Y(:,5) + pi/2 - atan2(Y(:,2),Y(:,1));
delta_vect = 5*parout.delta_vec;

n_start = 1;
step_size = 5;
n_end = 1500;

Xcg = 180;

% vector creator
k = 0;
TR_body_vect(1:length(n_start:step_size:n_end), 1) = {NaN(1,1)};
TR_engine_vect(1:length(n_start:step_size:n_end), 1) = {NaN(1,1)};
for ii = n_start:step_size:n_end
    k = k + 1;
    [TR_body, TR_engine] = animation_frame(theta_vect(ii), delta_vect(ii), Xcg);
    TR_body_vect{k} = TR_body;
    TR_engine_vect{k} = TR_engine;
end

% animation
figure, hold on, axis equal, grid on
for ii = 1:length(TR_engine_vect)
    fig1 = trisurf(TR_body_vect{ii}, 'FaceColor','black','EdgeColor','interp');
    fig2 = trisurf(TR_engine_vect{ii}, 'FaceColor','black','EdgeColor','interp');

    pause(0.05)
    drawnow

    view(0,0)

    delete(fig1)
    delete(fig2)
end

close all


%% Test

clear, clc
close all

theta = deg2rad(20);
delta = deg2rad(25);

[TR_body, TR_engine] = animation_frame(theta, delta);

figure
trisurf(TR_body, 'FaceColor','black','EdgeColor','interp');
axis equal, hold on, xlabel("x"), ylabel("y"), zlabel("z")
view(0,0)
trisurf(TR_engine, 'FaceColor','black','EdgeColor','interp');
