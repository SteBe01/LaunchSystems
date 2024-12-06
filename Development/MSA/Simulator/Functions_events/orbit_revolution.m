function [value, isterminal, direction] = orbit_revolution(t, y, params)
    persistent t_orbit
    if isempty(t_orbit)
        t_orbit = Inf;
    end

    T = 2*pi*sqrt((params.Re/1e3+400)^3/398600)*1.5;
    
    if norm(y(1:2)) >= 400e3+params.Re && t_orbit == Inf
        t_orbit = t;
    end

    value = t - (t_orbit+T);
    isterminal = 1;
    direction = 0;

end

