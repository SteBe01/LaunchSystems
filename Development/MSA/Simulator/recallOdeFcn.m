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
coeffs = zeros(length(T), 3);
dcm = zeros(2,2,length(T));
Thrust = zeros(length(T),1);

F_in = zeros(length(T), 2);
F_L_in = zeros(length(T), 2);
F_D_in = zeros(length(T), 2);

clear dyn
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
    coeffs(ii, :) = parout.coeffs;
    dcm(:,:,ii) = parout.dcm;

    F_in(ii, :) = parout.F_in;
    F_L_in(ii,:) = parout.F_L_in;
    F_D_in(ii,:) = parout.F_D_in;
    Thrust(ii,:) = parout.Thrust;
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
ode_out.coeffs = coeffs;
ode_out.dcm = dcm;

ode_out.F_in = F_in;
ode_out.F_L_in = F_L_in;
ode_out.F_D_in = F_D_in;
ode_out.Thrust = Thrust;

end

