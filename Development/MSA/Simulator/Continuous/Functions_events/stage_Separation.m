function [value, isterminal, direction] = stage_Separation(t, y, stage, params, varargin)
    persistent t_stop
    if isempty(t_stop)
        t_stop = Inf;
    end

    if ~isempty(varargin)
        [~,parout] = dyn(t, y, stage, params, 1, varargin{1});
    else
        [~,parout] = dyn(t, y, stage, params, 1);
    end
    if parout.F_in(2)/parout.m < 0 && parout.Thrust ~= 0
        exitFlag = 0;
    else
        exitFlag = 1;
    end

    if y(7) <= stage.m_prop_final && t_stop == Inf
        t_stop = t;
    end

    % value = y(2);
    value(1) = (t - (t_stop+3));
    value(2) = exitFlag;
    % value = t - (stage.t_burn_tot + stage.t_wait + 1);
    % value = y(2);
    isterminal = [1 1];
    direction = [0 0];
end