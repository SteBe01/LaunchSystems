function ode_out = recallOdeFcn(T, Y, stage, params, curr_stage)

qdyn = zeros(length(T), 1);
acc = zeros(length(T), 2);
alpha = zeros(length(T), 1);
moment = zeros(length(T), 1);
dv_drag_vec = zeros(length(T), 1);
rho = zeros(length(T), 1);
velsqq = zeros(length(T), 1);
m = zeros(length(T), 1);
dv_grav_vec = zeros(length(T), 1);
delta = zeros(length(T),1);

for ii = 1:length(T)
    [~, parout] = dyn(T(ii), Y(ii, :), stage, params, curr_stage);
    qdyn(ii) = parout.qdyn;
    acc(ii,:) = parout.acc;
    alpha(ii) = parout.alpha;
    moment(ii) = parout.moment;
    dv_drag_vec(ii) = parout.dv_drag;
    rho(ii) = parout.rho;
    velsqq(ii) = parout.velssqq;
    m(ii) = parout.m;
    dv_grav_vec(ii) = parout.dv_grav;
    delta(ii) = parout.delta;
end

ode_out.qdyn = qdyn;
ode_out.acc = acc;
ode_out.alpha = alpha;
ode_out.moment = moment;
ode_out.dv_drag_vec = dv_drag_vec;
ode_out.rho = rho;
ode_out.velsqq = velsqq;
ode_out.m = m;
ode_out.dv_grav_vec = dv_grav_vec;
ode_out.delta = delta;

end

