    function dxdt = fun_dyn(t, x, p, u)

% x(1) = theta
% x(2) = dottheta
% x(3) = alpha
% x(4) = zeta
% x(5) = dotzeta

dxdt = zeros(5,1);

dxdt(1) = x(2);
dxdt(2) = p.Ma*x(3) + p.Md*u;
dxdt(3) = -p.F/p.m/p.V*x(1) - p.N/p.m/p.V*x(3) + p.T/p.m/p.V*u + p.alphawdot;
dxdt(4) = x(5);
dxdt(5) = -p.F/p.m*x(1) - p.N/p.m*x(3) + p.T/p.m*u;

end