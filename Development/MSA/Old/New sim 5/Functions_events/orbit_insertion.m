function [value, isterminal, direction] = orbit_insertion(~, y)
    check = 1;
    if y(8) < 0
        check = 0;
    end

    value = y(4)*check;
    isterminal = 1;
    direction = 0;
end