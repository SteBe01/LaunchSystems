function ode_out = recallOdeFcn(T, Y, stage, params, curr_stage, delta)
    qdyn = zeros(length(T), 1);
    acc = zeros(length(T), 2);
    alpha = zeros(length(T), 1);
    moment = zeros(length(T), 1);
    for ii = 1:length(T)
        [~, parout] = dyn(T(ii), Y(ii, :), stage, params, curr_stage, delta);
        qdyn(ii) = parout.qdyn;
        acc(ii,:) = parout.acc;
        alpha(ii) = parout.alpha;
        moment(ii) = parout.moment;
    end

    ode_out.qdyn = qdyn;
    ode_out.acc = acc;
    ode_out.alpha = alpha;
    ode_out.moment = moment;
end