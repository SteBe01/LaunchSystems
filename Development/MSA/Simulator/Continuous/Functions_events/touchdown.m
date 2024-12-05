function [value, isterminal, direction] = touchdown(~, y, params)
    % value = y(2) - params.Re;
    value = 2;
    isterminal = 1;
    direction = 0;
end