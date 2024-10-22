clear
clc
close all

%% Preliminary launcher stages masses and volumes

% Data
m_pl  = 500;                  % Minimum to SSO, Requirement
DV    = 8165;                 % From Transfer.m, this one is from the Pegsus baseline to try and validate
g0    = 9.80665;
ve    = [290.2; 289.4; 287.4] * g0;      % Effective exhaust velocities, from typical Isp values (Pegasus Baseline)
eps_s = [0.0836; 0.0958; 0.141];        % Typical structural mass indexes  (Pegasus Baseline)



[m_vec, ~] = mass_compute(DV, m_pl, ve, eps_s);
GLOM       = sum(m_vec) + m_pl;     % Adding up (stage struct & prop) + payload
m_struct   = eps_s .* m_vec;        % By definition
m_prop     = m_vec - m_struct;      % By definition

avg_density = 1041.65;              % By data (Pegasus Baseline)
finessness  = 13.31;                % L/D=f (Pegasus Baseline)
V = GLOM / avg_density;             % V=m/rho
D = ( (4*V) / (pi*finessness) )^(1/3);  % V = L*pi*D^2/4 = pi*f*D^3/4   -> D^3 = 4*V/(pi*f)
L = finessness * D;


function [m_stage, lambda] = mass_compute(DV, m_pl, ve, eps_s)
% Optimizes stage masses using lagrangian multipliers method

% INPUTS:
% Target total Delta-v  [1x1]   [km/s]
% Payload mass          [1x1]   [kg]
% Eff. exhaust velocity [nx1]   [m/s]
% Structural mass index [nx1]   [-]

% Checks
if ~(isequal(size(DV),size(m_pl),[1,1])&&isequal(size(ve),size(eps_s)))
    error("Dimension of input data not correct")
end

% Initialization
options  = optimoptions('fsolve','Display','off');
n_stages = length(eps_s);
m_stage  = zeros(n_stages,1);

%Solving
f = @(lambda) DV - sum(ve.*log((ve*lambda-1)./(ve.*eps_s*lambda))); %Objective function

lambda = fsolve(f, 0.1, options);           % 0.1 is chosen at random as starting value different from zero
n = (ve*lambda-1)./(ve.*eps_s*lambda);      % The inverses of MR are computed

for i = n_stages:-1:1
    m_stage(i) = (n(i)-1)/(1-n(i)*eps_s(i)) * (m_pl+sum(m_stage));   %This mass is structure and propellant of indicated stage
end

end