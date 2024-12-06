function [dY] = dyn(t,Y,stage,params)

% Data extraction
x = Y(1);
y = Y(2);
altitude = -Y(3);
x_dot = Y(4);
y_dot = Y(5);
z_dot = Y(6);
p = Y(7);
q = Y(8);
r = Y(9);
qw = Y(10);
qx = Y(11);
qy = Y(12);
qz = Y(13);

vels_body = [x_dot y_dot z_dot]';
velsNorm = norm(vels_body);

Q = [qw qx qy qz];
Q = Q/norm(Q);
dcm = quat2dcm(Q);

vels_intertial = dcm'*vels_body;



t_wait = 0;

% Retrieve data used multiple times 
t_burn_tot = stage.t_burn_tot;
Re = params.Re;

% Mass estimation
if t <= 1000
    m = stage.m0;
elseif t > t_wait && t <= t_burn_tot + t_wait
    m = stage.m0 - stage.m_dot * (t-t_wait);
else
    m = stage.m0 - stage.m_dot * t_burn_tot;
end

% Forces
rho = getDensity(altitude);
qdyn = 0.5*rho*velsNorm^2;
g = params.g0/((1+altitude/Re)^2);
S = pi*(stage.d^2/4);

if t<0
    T = stage.Thrust;
else
    T = 0;
end

CA = stage.Cd;
CY = 0;
CN = stage.Cl;

fX = qdyn*S*CA;
fY = qdyn*S*CY;
fZ = qdyn*S*CN;
Fg = dcm*[0; 0; m*g];
F = Fg + [-fX+T, fY, -fZ]';




du = F(1)/m - q*z_dot + r*y_dot;
dv = F(2)/m - r*x_dot + p*z_dot;
dw = F(3)/m - p*y_dot + q*x_dot;

% Rotation
% dp = (Iyy - Izz)/Ixx*q*r + qdynL_V/Ixx*(velsNorm*Cl+Clp*p*C/2) - Ixxdot*p/Ixx;
% dq = (Izz - Ixx)/Iyy*p*r + qdynL_V/Iyy*(velsNorm*Cm + (Cmad+Cmq)*q*C/2)...
%     - Iyydot*q/Iyy;
% dr = (Ixx - Iyy)/Izz*p*q + qdynL_V/Izz*(velsNorm*Cn + (Cnr*r+Cnp*p)*C/2)...
%     - Izzdot*r/Izz;

Ixx = stage.I;
Iyy = stage.J;
Izz = stage.K;

Ixxdot = 0;
Iyydot = 0;
Izzdot = 0;

Cl = 0.5;
Clp = 0.5;
Cmq = 0.5;
Cmad = 0.5;
Cm = 0.5;
Cnr = 0.5;
Cn = 0.5;
Cnp = 0.5;

C = stage.d;
qdynL_V = 0.5*rho*velsNorm*S*C;

dp = (Iyy - Izz)/Ixx*q*r + qdynL_V/Ixx*(velsNorm*Cl+Clp*p*C/2) - Ixxdot*p/Ixx;
dq = (Izz - Ixx)/Iyy*p*r + qdynL_V/Iyy*(velsNorm*Cm + (Cmad+Cmq)*q*C/2) - Iyydot*q/Iyy;
dr = (Ixx - Iyy)/Izz*p*q + qdynL_V/Izz*(velsNorm*Cn + (Cnr*r+Cnp*p)*C/2) - Izzdot*r/Izz;



% Quaternions
OM = [ 0 -p -q -r
    p  0  r -q
    q -r  0  p
    r  q -p  0];
dQQ = 1/2*OM*Q';

% Derivatives
dY = zeros(13, 1);

dY(1:3) = vels_intertial;
dY(4:6) = [du; dv; dw];
dY(7:9) = [dp; dq; dr];
dY(10:13) = dQQ;

end

