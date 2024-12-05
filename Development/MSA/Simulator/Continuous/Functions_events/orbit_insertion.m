function [value, isterminal, direction] = orbit_insertion(~, y)
    % value = 1;
    % value = norm(y(1:2)) - (400e3+6378e3);
    value = y(4);
    isterminal = 1;
    direction = 0;
end