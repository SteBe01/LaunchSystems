function [TR_body, TR_engine] = animation_frame(theta, delta, Xcg)

TR_body = stlImport_internal("Parte1.STL");
TR_engine = stlImport_internal("Parte2.STL");

angle1 = deg2rad(0);
angle2 = deg2rad(0);
angle3 = deg2rad(90);

R1 = @(ang) [1 0 0
            0 cos(ang) -sin(ang)
            0 sin(ang) cos(ang)];
R2 = @(ang) [cos(ang) 0 sin(ang)
            0 1 0
            -sin(ang) 0 cos(ang)];
R3 = @(ang) [cos(ang) -sin(ang) 0
            sin(ang) cos(ang) 0
            0 0 1];
R_tot = @(angle1, angle2, angle3) R3(angle1)*R2(angle2)*R1(angle3);

A = R_tot(angle1, angle2, angle3);
[TR_body] = stlRotate_internal(TR_body, A);
[TR_engine] = stlRotate_internal(TR_engine, A);

TR_body = stlShift(TR_body, max(TR_body.Points(:,1)));
TR_engine = stlShift(TR_engine, max(TR_engine.Points(:,1)));

% Engine rotation
A = R_tot(0, -delta, 0);
[TR_engine] = stlRotate_internal(TR_engine, A);

% Engine shift
TR_engine = stlShift(TR_engine, max(TR_body.Points(:,1))-min(TR_body.Points(:,1)));

% Xcg shift
TR_body = stlShift(TR_body, -Xcg);
TR_engine = stlShift(TR_engine, -Xcg);

% Whole attitude rotation
A = R_tot(0, -theta, 0);
[TR_body] = stlRotate_internal(TR_body, A);
[TR_engine] = stlRotate_internal(TR_engine, A);


%% Functions

function [TR] = stlImport_internal(stl_file_name)
    gm = stlread(stl_file_name);
    
    x_dim = max(gm.Points(:,1))-min(gm.Points(:,1));
    y_dim = max(gm.Points(:,2))-min(gm.Points(:,2));
    z_dim = max(gm.Points(:,3))-min(gm.Points(:,3));
    
    % centre of the body
    cm_body = [x_dim/2 y_dim/2 z_dim/2];
    transposition = [min(gm.Points(:,1)) min(gm.Points(:,2)) min(gm.Points(:,3))];
    
    % points with centre of the body in [0 0 0]
    new_points = gm.Points - ones(size(gm.Points)).*(cm_body + transposition);
    TR = triangulation(gm.ConnectivityList, new_points);
end

function [TR] = stlShift(TR, shift)
    transposition = [shift 0 0];
    
    % points with centre of the body in [0 0 0]
    new_points = TR.Points - ones(size(TR.Points)).*(transposition);
    TR = triangulation(TR.ConnectivityList, new_points);
end

function [TR] = stlRotate_internal(TR, Rot_matrix)
    newTrPoints = zeros(size(TR.Points));
    for i = 1:size(TR.Points, 1)
        newTrPoints(i,:) = Rot_matrix * TR.Points(i,:)';
    end
    TR = triangulation(TR.ConnectivityList, newTrPoints);
end


end

