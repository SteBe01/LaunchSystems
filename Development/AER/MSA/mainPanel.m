% clear
% close all
% clc
%% MAIN

R_v = [         % vettore con i vari raggi (ogiva, corpo, spalla)
    0           % raggio iniziale = 0m (punta)
    0.6         % raggio dell'ogiva (finale)
    0.6         % raggio della spalla (finale)
    0.6         % raggio finale flare/boatt
    ];

L_v =[
    2.4           % lunghezza ogiva
    4.2876          % lunghezza corpo
    4.2876           % lunghezza spalla/boatt
    4.2876          % lunghezza corpo
    4.2876           % lunghezza flare/boatt
    ];



geometrySegments = struct(...
    'xStart', {0, L_v(1), sum(L_v(1:2)), sum(L_v(1:3)), sum(L_v(1:4))}, ...
    'xEnd', {L_v(1), sum(L_v(1:2)), sum(L_v(1:3)), sum(L_v(1:4)), sum(L_v(1:5))},...
    'funcZ', {@(x) R_v(2)/L_v(1) * x, ...            % --> NON SERVE: mettere in input solo funcR
    @(x) R_v(2) + 0 * x, ...
    @(x) R_v(2) + (R_v(3)-R_v(2))/L_v(3) * (x - sum(L_v(1:2))), ...
    @(x) R_v(3) + 0 * x, ...
    @(x) R_v(3) + (R_v(4)-R_v(3))/L_v(5) * (x - sum(L_v(1:4)))}, ...
    'funcR', {@(x) -R_v(2)/L_v(1)^2 .* x.^2 + 2 * R_v(2)/L_v(1) .* x, ...
    @(x) R_v(2) + 0 * x, ...
    @(x) R_v(2) + (R_v(3)-R_v(2))/L_v(3) * (x - sum(L_v(1:2))), ...
    @(x) R_v(3) + 0 * x, ...
    @(x) R_v(3) + (R_v(4)-R_v(3))/L_v(5) * (x - sum(L_v(1:4)))});
%%%%%%%%%%%%%%%% funcR for TANGENT OGIVE %%%%%%%%%%%%%%%%%%
    % 'funcR', {@(x) -R_v(2)/L_v(1)^2 .* x.^2 + 2 * R_v(2)/L_v(1) .* x, ...
    % @(x) R_v(2) + 0 * x, ...
    % @(x) R_v(2) + (R_v(3)-R_v(2))/L_v(3) * (x - sum(L_v(1:2))), ...
    % @(x) R_v(3) + 0 * x, ...
    % @(x) R_v(3) + (R_v(4)-R_v(3))/L_v(5) * (x - sum(L_v(1:4)))});

    % 'funcR', {@(x) R_v(2)/L_v(1) * x, ...                   % CONICAL
    % @(x) R_v(2) + 0 * x, ...
    % @(x) R_v(2) + (R_v(3)-R_v(2))/L_v(3) * (x - sum(L_v(1:2))), ...
    % @(x) R_v(3) + 0 * x, ...
    % @(x) R_v(3) + (R_v(4)-R_v(3))/L_v(5) * (x - sum(L_v(1:4)))});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


D_ref = 2 * max(R_v);
S_ref = D_ref^2 * pi/4;     % misure di riferimento per adimensionalizzazione

% alpha = 10;
% Mach = 6;


altitude = 20000;
numPanels_LONG = 10;        % for each section of the body
numPanels_RAD = 30;         % each long station ther will be discretised by # points


% [CL_New, CD_New, CL_ModNew, CD_ModNew] = panelMethodRocket3D(geometrySegments, alpha, Mach, altitude, numPanels_LONG, numPanels_RAD, L_v, S_ref)



% CICLO SU ALPHA E MACH:

alpha_v = 0:5:10;      % in gradi!!
Mach_v = 3:8;

CL_New = [];
CD_New = [];
CL_ModNew = [];
CD_ModNew = [];

for j = 1:length(Mach_v)
    for i = 1:length(alpha_v)
    
        if alpha_v(i) > 90 || alpha_v(i) < 0
            fprintf('alpha maggiore di 90 gradi!!, deve essere compreso tra 0 e 90')
            break;
        end
    
        close all
        [CL_New(j, i), CD_New(j, i), CL_ModNew(j, i), CD_ModNew(j, i)] = panelMethodRocket3D(geometrySegments, alpha_v(i), Mach_v(j), altitude, numPanels_LONG, numPanels_RAD, L_v, S_ref);
    
    end
end

figure
plot(Mach_v, CL_New)
title('CL New vs Mach')
figure
plot(Mach_v, CD_New)
title('CD New vs Mach')
figure
plot(Mach_v, CL_ModNew)
title('CL ModNew vs Mach')
figure
plot(Mach_v, CD_ModNew)
title('CD ModNew vs Mach')


save('CL_NEW.mat', 'CL_New')
save('CL_MODNEW.mat', 'CL_ModNew')
save('CD_NEW.mat', 'CD_New')
save('CD_MODNEW.mat', 'CD_ModNew')






%%
function [CL_New, CD_New, CL_ModNew, CD_ModNew] = panelMethodRocket3D(geometrySegments, alpha, Mach, altitude, numPanels_LONG, numPanels_RAD, L_v, S_ref)
    % INPUTS:
    % geometrySegments: Struct array defining segments of the rocket geometry
    % alpha: Angle of attack in degrees
    % Mach: Freestream Mach number
    % altitude: Altitude in meters
    % numPanels_LONG: Number of panels for the surface discretization longitudinally
    % numPanels_RAD: # panels for discretiz at fixed long station

    % Constants
    gamma = 1.4; % Ratio of specific heats for air

    % Atmospheric model (ISA standard atmosphere)
    [T, a, P, rho] = atmosisa(altitude);

    % Freestream properties
    V = Mach * a; % Freestream velocity
    q = 0.5 * rho * V^2;    % Pressione dinamica
    P01 = P * (1 + (gamma - 1) / 2 * Mach^2)^(gamma / (gamma - 1));      % total preassure free stream
    P02 = P01 * ((gamma + 1) * Mach^2 / ((gamma - 1) * Mach^2 + 2))^(gamma / (gamma - 1)) * ...
        ( (gamma + 1) / ( 2 * gamma * Mach^2 - (gamma - 1) ) )^(1/(gamma - 1));

    % Discretize the complex geometry
    [x, z, r] = discretizeAxisymmetricGeometry(geometrySegments, numPanels_LONG, L_v);

    
    % Visualize the 3D discretized geometry
    [X, Y, Z, panelNormals, panelAreas] = GeometryWithPanels(x, z, r, numPanels_RAD);

    
    % NEWTON METHODS:
    % Convert alpha to radians
    alphaRad = deg2rad(alpha);


    % Calculate aero Force using Newton's method:
    [F_aero_New] = NewtonSolver(alphaRad, q, panelNormals, panelAreas);
    F_assiale_New = sum(F_aero_New(:, 1))
    F_normale_New = sum(F_aero_New(:, 3))

    % Passo a coefficienti e assi vento:
    Cn_New = F_normale_New / (q * S_ref);
    Ca_New = F_assiale_New / (q * S_ref);
    [CL_New, CD_New] = body2wind(Cn_New, Ca_New, alpha);



    % Calculate aero Force using MODIFIED Newton's method:
    [F_aero_ModNew] = ModNewtonSolver(alphaRad, q, panelNormals, panelAreas, P02, P);
    F_assiale_ModNew = sum(F_aero_ModNew(:, 1))
    F_normale_ModNew = sum(F_aero_ModNew(:, 3))

    % Passo a coefficienti e assi vento:
    Cn_ModNew = F_normale_ModNew / (q * S_ref);
    Ca_ModNew = F_assiale_ModNew / (q * S_ref);
    [CL_ModNew, CD_ModNew] = body2wind(Cn_ModNew, Ca_ModNew, alpha);

end



%%
function [x, z, r] = discretizeAxisymmetricGeometry(geometrySegments, numPanels_LONG, L_v)
    % Discretizes the geometry into panels for axisymmetric calculations
    x = [];
    z = [];
    r = [];
    for segment = geometrySegments

        if segment.xEnd == sum(L_v)
            segmentX = linspace(segment.xStart, segment.xEnd, numPanels_LONG);
        else
            segmentX = linspace(segment.xStart, segment.xEnd, numPanels_LONG);
            segmentX = segmentX(1:(length(segmentX)-1));
        end

        segmentZ = segment.funcZ(segmentX);
        segmentR = segment.funcR(segmentX);
        x = [x, segmentX];
        z = [z, segmentZ];
        r = [r, segmentR];
    end

    % z = r;

end



%%
function [X, Y, Z, panelNormals, panelAreas] = GeometryWithPanels(x, z, r, numPanels_RAD)
    % Creates a 3D visualization of the axisymmetric geometry with panels
    theta = linspace(0, 2 * pi, numPanels_RAD); % Azimuthal angle for rotation
    X = repmat(x', 1, length(theta));
    Y = -r' .* cos(theta);
    Z = -r' .* sin(theta);

    % Visualize geometry
    figure;
    surf(X, Y, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    colormap turbo;
    % shading interp;
    hold on;
    axis equal

    % Calculate and visualize panel normals
    [panelNormals, panelAreas, centerX, centerY, centerZ] = calculatePanelNormals(X, Y, Z);

    % Visualize panel vertices
    scatter3(X(:), Y(:), Z(:), 'filled', 'MarkerEdgeColor', 'k', ...
             'MarkerFaceColor', 'b', 'DisplayName', 'Panel Vertices', 'LineWidth', 0.1);

    % Visualize panel normals
    quiver3(centerX, centerY, centerZ, ...
            squeeze(panelNormals(:, :, 1))/2, ...
            squeeze(panelNormals(:, :, 2))/2, ...
            squeeze(panelNormals(:, :, 3))/2, ...
            'k', 'LineWidth', 1, 'DisplayName', 'Panel Normals', 'AutoScale','off');

    % Plot settings
    title('3D Geometry with Panel Normals and Vertices');
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    axis equal;
    grid on;
    legend show; % Add legend to distinguish vertices and normals
end



%% NEWTON:
function [F_aero] = NewtonSolver(alphaRad, q, panelNormals, panelAreas)
    
    F_aero = [];

    % !! ALPHA in RADIANS !!

    for i = 2:size(panelNormals, 1)     % scorre le stazioni longitudinali
        for j = 2:size(panelNormals, 2) % scorre le theta per ogni stazione

            normal = [panelNormals(i, j, 1), panelNormals(i, j, 2), panelNormals(i, j, 3)];
            beta = angleBetweenPanelAndFlow(normal, alphaRad);

            Cp = 2 * (sin(beta))^2;
            F_aero = [F_aero;  - Cp * q * panelAreas(i, j) * normal];
        end
    end

end



%% MODIFIED NEWTON:
function [F_aero_ModNew] = ModNewtonSolver(alphaRad, q, panelNormals, panelAreas, P02, P)
    
    F_aero_ModNew = [];

    % !! ALPHA in RADIANS !!

    for i = 2:size(panelNormals, 1)     % scorre le stazioni longitudinali
        for j = 2:size(panelNormals, 2) % scorre le theta per ogni stazione

            normal = [panelNormals(i, j, 1), panelNormals(i, j, 2), panelNormals(i, j, 3)];
            beta = angleBetweenPanelAndFlow(normal, alphaRad);

            Cp_max = (P02 - P) / (q);
            Cp = Cp_max * (sin(beta))^2;
            F_aero_ModNew = [F_aero_ModNew;  - Cp * q * panelAreas(i, j) * normal];
        end
    end

end



%%
function [CL, CD] = calculateAerodynamicCoefficients(Cp, lengths, radii, nx, nz, alphaRad)
    % Calculate aerodynamic coefficients CL and CD for axisymmetric body
    S_ref = 2 * pi * trapz(radii .* lengths); % Reference area
    Fx = -Cp .* lengths .* radii .* nx; % Axial force contribution
    Fz = -Cp .* lengths .* radii .* nz; % Normal force contribution
    CL = sum(Fz * cos(alphaRad) - Fx * sin(alphaRad)) / S_ref;
    CD = sum(Fx * cos(alphaRad) + Fz * sin(alphaRad)) / S_ref;
end



%%
function [panelNormals, panelAreas, centerX, centerY, centerZ] = calculatePanelNormals(X, Y, Z)
    % INPUT:
    % X, Y, Z: Coordinates of panel vertices from GeometryWithPanels()
    % OUTPUT:
    % panelNormals: Normal vectors of each panel
    % centerX, centerY, centerZ: Centers of the panels for visualization

    numPanelsX = size(X, 1) - 1; % Number of panels along X
    numPanelsTheta = size(X, 2) - 1; % Number of panels along theta

    panelNormals = zeros(numPanelsX, numPanelsTheta, 3); % Initialize normals
    panelAreas = zeros(numPanelsX, numPanelsTheta);
    centerX = zeros(numPanelsX, numPanelsTheta);
    centerY = zeros(numPanelsX, numPanelsTheta);
    centerZ = zeros(numPanelsX, numPanelsTheta);

    % Loop through each panel
    for i = 1:numPanelsX
        for j = 1:numPanelsTheta
            % Panel vertices
            v1 = [X(i, j), Y(i, j), Z(i, j)];
            v2 = [X(i + 1, j), Y(i + 1, j), Z(i + 1, j)];
            v3 = [X(i, j + 1), Y(i, j + 1), Z(i, j + 1)];
            v4 = [X(i + 1, j + 1), Y(i + 1, j + 1), Z(i + 1, j + 1)];

            % figure
            % scatter3(v1(1), v1(2), v1(3))
            % hold on
            % scatter3(v2(1), v2(2), v2(3))
            % scatter3(v3(1), v3(2), v3(3))
            % scatter3(v4(1), v4(2), v4(3))
            % legend
            % Edge vectors
            edge1 = v2 - v3;
            edge2 = v4 - v1;

            % Compute normal as cross product of edges
            normal = cross(edge1, edge2);
            panelNormals(i, j, :) = - normal / norm(normal); % Normalize

            % Compute panel Area:
            panelAreas(i, j) = norm(normal) * 0.5;

            % Compute panel center
            centerX(i, j) = mean([v1(1), v2(1), v3(1), v4(1)]);
            centerY(i, j) = mean([v1(2), v2(2), v3(2), v4(2)]);
            centerZ(i, j) = mean([v1(3), v2(3), v3(3), v4(3)]);

        end
    end
end



%% ANGOLO ATTACCO del singolo pannello:
function beta = angleBetweenPanelAndFlow(normal, alpha)
    % normal: [nx, ny, nz], unit normal vector of the panel
    % alpha: flow angle in the xz plane (in RADIANS)

    % Flow vector in the xz plane
    flow = [cos(alpha), 0, sin(alpha)]; % Unit flow vector

    % Dot product between flow and normal
    dotProduct = dot(normal, flow);

    % beta = acos(dotProduct/(norm(flow) * norm(normal)));
    % beta = beta - pi/2;

    if dotProduct > 0
        beta = 0;
    else
        % Angle between normal and flow
        % beta = 90 - dotProduct;
        beta = acos(dotProduct/(norm(flow) * norm(normal))); % In radians
        beta = beta - pi/2;
    end
end



%% Body to Wind:
function [CL, CD] = body2wind(Cn, Ca, alpha_v)

CL = [];
CD = [];

CL = Cn .* cos(deg2rad(alpha_v)) - Ca .* sin(deg2rad(alpha_v));
CD = Cn .* sin(deg2rad(alpha_v)) + Ca .* cos(deg2rad(alpha_v));
end