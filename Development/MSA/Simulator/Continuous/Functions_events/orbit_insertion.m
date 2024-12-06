function [value, isterminal, direction] = orbit_insertion(~, y)
    % value = 1;
    % value = norm(y(1:2)) - (400e3+6378e3);

    beta = atan2(y(2), y(1));
    ang = pi/2-beta;
    vec_rotated = [cos(ang) -sin(ang); sin(ang) cos(ang)]*y(3:4);

    % value = vec_rotated(2);
    value = 1;
    isterminal = 1;
    direction = 0;
end