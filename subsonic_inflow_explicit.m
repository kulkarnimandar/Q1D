function prim_cc_gi = subsonic_inflow_explicit(q1,T0,p0,CSA_flag,prim_cc_gi_flow)

if nargin == 3
    CSA_flag = 0;
end

prim_cc_gi = zeros(3,1);
g = 1.4;
R = 287;

if CSA_flag ~= 1
    u_gi = q1(2);
    prim_cc_gi(2) = u_gi;
    
    T = T0 - (g-1)/(2*g*R)*u_gi^2;
    psi = T0/T;
    p = p0/psi^(g/(g-1));
    rho0 = p0/(R*T0);
    rho = rho0/psi^(1/(g-1));
    prim_cc_gi(1) = rho;
    prim_cc_gi(3) = p;
else
    % For Convergent only nozzle with exit at x = -0.1
%     prim_cc_gi =   [ 0.029681158903921;
%                     -85.967594122473116
%                      7141.905060096178];

% For C-D nozzle
% prim_cc_gi = [0.048278639525954;
%               -1.172102369895277e+002;
%               1.160738552859373e+004];
%           prim_cc_gi = 2*q1 - q2;

    u_gi = prim_cc_gi_flow(2);
    T_gi = T0 - (g-1)/(2*g*R)*u_gi^2;
    psi_gi = T0/T_gi;
    
    
    u_gi_prime = q1(2);
    prim_cc_gi(2) = u_gi_prime;
    
    T_gi_prime = - (g-1)/(2*g*R)*2*u_gi*u_gi_prime;
    psi_gi_prime = -T0/T_gi^2*T_gi_prime;
    p_gi_prime = (-g/(g-1)) * p0/psi_gi^(1 + g/(g-1)) * psi_gi_prime;
    rho0 = p0/(R*T0);
    rho_gi_prime = (-1/(g-1)) * rho0/psi_gi^(1 + 1/(g-1)) * psi_gi_prime;
    
    prim_cc_gi(1) = rho_gi_prime;
    prim_cc_gi(3) = p_gi_prime;
end
end