function [value, isterminal, direction] = touchdown(~, y, params)
    value = norm(y(1:2)) - params.Re;
    isterminal = 1;
    direction = 0;
end