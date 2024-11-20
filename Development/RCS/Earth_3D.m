function Earth_3D(Re,alpha,r_n)

% Earth_3D.m - Earth texture loaded in a plot.
%
% PROTOTYPE:
%   Earth_3D(Re,alpha,r_n)
%
% DESCRIPTION:
%   Function to load the Earth modelled as a sphere inside a figure.
%
% INPUT:
%   Re          [1x1]       Earth mean radius       [km]
%   alpha       [1x]        Set the transparency of the globe: 1 = opaque, 0 = invisible 
%   r_n         [3x1]       Position vector wrt Sun [km]
%
% OUTPUT:
%   []          [figure]    Figure open with the Earth picture loaded
%
% CONTRIBUTORS:
% Camilla Airoldi, Jesus Braza Rodriguez, Elena De Marco, Riccardo Monti
%
% VERSIONS
% 2023-12-01 Last version

hold on;

% Set the axes scale equal
axis equal;

% Set initial view
view(120,30);

% Define the number of panels to be used to model the sphere 
npanels = 180;  

if nargin == 1
[x, y, z] = ellipsoid(0, 0, 0, Re, Re, Re, npanels);
else 
[x, y, z] = ellipsoid(r_n(1,1), r_n(2,1), r_n(3,1), Re, Re, Re, npanels);
end

if nargin == 1
    alpha = 1;
end

% Create the globe with the surf function
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 'none');

% Load Earth image for texture map
cdata = imread("Earth_image.jpg");

% Set the transparency of the globe: 1 = opaque, 0 = invisible 

set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');

end