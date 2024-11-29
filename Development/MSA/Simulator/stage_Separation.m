function [value, isterminal, direction] = stage_Separation(t, y, stage)
    persistent t_stop
    if isempty(t_stop)
        t_stop = Inf;
    end

    if y(7) <= stage.m_prop_final && t_stop == Inf
        t_stop = t;
    end

    value = t - (t_stop+3);
    % value = t - (stage.t_burn_tot + stage.t_wait + 1);
    % value = y(2);
    isterminal = 1;
    direction = 0;
end