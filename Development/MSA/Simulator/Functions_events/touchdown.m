function [value, isterminal, direction] = touchdown(~, y, params)
    value = y(2) - params.Re;
    isterminal = 1;
    direction = 0;
end