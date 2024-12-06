function [value, isterminal, direction] = stop_ev(~, y)
    value = y(8);
    isterminal = 1;
    direction = 0;
end

