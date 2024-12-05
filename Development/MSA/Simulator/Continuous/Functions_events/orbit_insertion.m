function [value, isterminal, direction] = orbit_insertion(~, y)
    value = y(4);
    isterminal = 1;
    direction = 0;
end