function [L] = L(lambda, r, s, f, h)

% lambda = taper ratio
% r = radius of body at tail[m]
% s = maximum semispsan of tail in combination with body [m]
% f = wing vortex semispan at wing trailing edge [m]
% h = height of wing vortex above body axis at tail centre of pressure [m]


L = (((s-lambda*r)-f*(1-lambda))/(2*(s-r))*log((h^2+(f-s)^2)/(h^2+(f-r)^2))-...
    (1-lambda)/(s-r)*((s-r)+h*atan((f-s)/h)-h*atan((f-r)/h)));



end







