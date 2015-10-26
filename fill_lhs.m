function [L, D, U, DGi, UGi, LGo, DGo] = fill_lhs(flux_scheme, dx, area_f, dadx_cc, prim_cc, prim_gc, pback)

cells = length(dx);
faces = cells+1;
L = zeros(3,3,cells);
D = zeros(3,3,cells);
U = zeros(3,3,cells);
prim_cc_gi = prim_gc(:,1);
prim_cc_go = prim_gc(:,2);
jac_function = str2func(['jac_',flux_scheme]);
%% Fill L, D, U from jac_L and jac_R
for f = 2:faces-1
    c = f-1; % so c=1:faces-2, means c=1:cells-1, because faces=cells+1
    [jac_L,jac_R] = jac_function(prim_cc(:,c),prim_cc(:,c+1));
    % (a)
    D(:,:,c) = D(:,:,c) + jac_L*area_f(f);
    U(:,:,c) = U(:,:,c) + jac_R*area_f(f);
    % Contribution of f1 to L(:,:,1) and D(:,:,1) is missing
    % (b)
    L(:,:,c+1) = L(:,:,c+1) - jac_L*area_f(f);
    D(:,:,c+1) = D(:,:,c+1) - jac_R*area_f(f);
    % Contribution of f_(N+1) to D(:,:,N) and U(:,:,N) is missing
end
% Contribution of f1 to D(:,:,1) and L(:,:,1)
[jac_L,jac_R] = jac_function(prim_cc_gi,prim_cc(:,1));
L(:,:,1) = L(:,:,1) - jac_L*area_f(1);
D(:,:,1) = D(:,:,1) - jac_R*area_f(1);
% Contribution of f_(N+1) to D(:,:,N) and U(:,:,N)
[jac_L,jac_R] = jac_function(prim_cc(:,cells),prim_cc_go);
D(:,:,cells) = D(:,:,c) + jac_L*area_f(faces);
U(:,:,cells) = U(:,:,c) + jac_R*area_f(faces);

%% Add source Jacobian term
for c = 1:cells
    source_jac = jac_source_1D(dadx_cc(c));
    D(:,:,c) = D(:,:,c) - source_jac*dx(c);
end

%% DGi and UGi for the input ghost cell
% DGi, UGi are directly on the top of L1, D1
% Res_gi = [psi_1;
%           psi_2;
%           psi_3]; 
%
% DGi = [dpsi_1_drho0, dpsi_1_du0, dpsi_1_dp0;
%        dpsi_2_drho0, dpsi_2_du0, dpsi_2_dp0;
%        dpsi_3_drho0, dpsi_3_du0, dpsi_3_dp0];
%    
% UGi = [0,      0,     0;
%        0, dpsi_2_du1, 0;
%        0,      0,     0];
g = 1.4;
R = 287;
%
rho0 = prim_cc_gi(1);
u0 = prim_cc_gi(2);
p0 = prim_cc_gi(3);
% u1 = prim_cc(2,1);
%
T0 = p0/(rho0*R);
a0 = speed_of_sound(p0, rho0);
M0 = u0/a0;
phi = 1 + (g-1)/2*M0^2;
%
dT0_drho0 = -p0/(rho0^2*R);
dT0_du0 = 0;
dT0_dp0 = 1/(rho0*R);
dphi_drho0 = (g-1)/(2*g)*u0^2/p0;
dphi_du0 = (g-1)/(2*g)*rho0*2*u0/p0;
dphi_dp0 = (g-1)/(2*g)*rho0*u0^2*(-1/p0^2);
% T_stag = T0*phi;
% psi_1 = T_stag - T_stag_BC;
dT_stag_drho0 = T0*dphi_drho0 + phi*dT0_drho0;
dT_stag_du0 = T0*dphi_du0 + phi*dT0_du0;
dT_stag_dp0 = T0*dphi_dp0 + phi*dT0_dp0;
dpsi_1_drho0 = dT_stag_drho0;
dpsi_1_du0 = dT_stag_du0;
dpsi_1_dp0 = dT_stag_dp0;
%
% psi_2 = u0 - u1;
dpsi_2_drho0 = 0;
dpsi_2_du0 = 1;
dpsi_2_dp0 = 0;
dpsi_2_du1 = -1;
%
% P_stag = p0*phi^(g/(g-1));
% psi_3 = P_stag - P_stag_BC;
dP_stag_drho0 = p0*(g/(g-1))*phi^(g/(g-1) - 1)*dphi_drho0;
dP_stag_du0   = p0*(g/(g-1))*phi^(g/(g-1) - 1)*dphi_du0;
dP_stag_dp0   = p0*(g/(g-1))*phi^(g/(g-1) - 1)*dphi_dp0 + phi^(g/(g-1));
dpsi_3_drho0 = dP_stag_drho0;
dpsi_3_du0 = dP_stag_du0;
dpsi_3_dp0 = dP_stag_dp0;
%
DGi = [dpsi_1_drho0, dpsi_1_du0, dpsi_1_dp0;
       dpsi_2_drho0, dpsi_2_du0, dpsi_2_dp0;
       dpsi_3_drho0, dpsi_3_du0, dpsi_3_dp0];
%
UGi = [0,      0,     0;
       0, dpsi_2_du1, 0;
       0,      0,     0];

%% LGo, DGo for the output ghost cell
% LGo, DGo are directly on below DN, UN
if pback > 0
    % If pback specified, for subsonic outflow case, p_(N+1) = p_back,
    % Hence dpsi_3_dpN = 0, So, last line in the matrices changes.
    % rho and u are still extrapolated
    LGo = -1*eye(3,3);
    DGo = eye(3,3);
    LGo(3,3) = 0;
else
    % Supersonic outflow BC is based on extrapolating from interior cells
    LGo = -1*eye(3,3);
    DGo = eye(3,3);
end

end

function [jac_l,jac_r] = jac_flux_vanLeer(qL,qR)
g =1.4;
% handle left side first

% store off primitive variables from vector
    rho = qL(1);
    u = qL(2);
    p = qL(3);
    rhoet = p*(1/(g-1)) + (1/2) * rho * u^2;
    a = speed_of_sound(p, rho);
    m = u/a;

% linearization of primitive variables wrt primitive variables
    drho_dq = [1 0 0];

    du_dq = [0 1 0];

    dp_dq = [0 0 1];

    drhoet_dq = [(1/2)*u^2, rho*u, (1/(g-1))];

% speed of sound linearization
    da_dq = [ -(1/2)*a/rho, 0, (1/2)*g/(rho*a) ];

% mach number linearization
    dm_dq = du_dq/a - (u/a^2)*da_dq;

% linearization of FVS mass term
    fa = (1/4)*rho*a*(m+1)^2;

    dfa_dq =  (1/4)*( drho_dq*a*(m+1)^2 + ...
                 rho*da_dq*(m+1)^2 + rho*a*2*(m+1)*dm_dq );

% linearization of FVS energy term
    fb = a*((g-1)*m + 2);

    dfb_dq = da_dq*((g-1)*m + 2) + dm_dq*(g-1)*a;

% Flux = [ fa, fa*fb*xg, (1/2)*fa*fb*fb*xg2m1 ]

% now, formulate jacobian

    jac_l = zeros(3,3);

    if ( abs(m) < 1 )
      jac_l(1,:) = dfa_dq;
      jac_l(2,:) = (dfa_dq*fb + dfb_dq*fa)*(1/g);
      jac_l(3,:) = (1/2)*(1/(g^2-1))*(dfa_dq*fb*fb + 2*fa*fb*dfb_dq);
    elseif ( m >= 1 )
      jac_l(1,:) = drho_dq*u + rho*du_dq;
      jac_l(2,:) = drho_dq*u*u + 2*rho*u*du_dq + dp_dq;
      jac_l(3,:) = ( drhoet_dq + dp_dq )*u + (rhoet + p)*du_dq;
    end

% handle right side second

% store off primitive variables from vector
    rho = qR(1);
    u = qR(2);
    p = qR(3);
    rhoet = p*(1/(g-1)) + (1/2) * rho * u^2;
    a = speed_of_sound(p,rho);
    m = u/a;

% linearization of right primitive variables wrt primitive variables
    drho_dq = [1 0 0];

    du_dq = [0 1 0];

    dp_dq = [0 0 1];

    drhoet_dq = [(1/2)*u^2, rho*u, (1/(g-1))];

% speed of sound linearization
    da_dq = [-(1/2)*a/rho, 0, (1/2)*g/(rho*a)];

% mach number linearization
    dm_dq = du_dq/a - (u/a^2)*da_dq;

% linearization of FVS mass term
    fa = -(1/4)*rho*a*(m-1)^2;

    dfa_dq = -(1/4)*( drho_dq*a*(m-1)^2 + ...
             rho*da_dq*(m-1)^2 + rho*a*2*(m-1)*dm_dq );

% linearization of FVS energy term

    fb = a*((g-1)*m - 2);

    dfb_dq = ((g-1)*m - 2)*da_dq + (g-1)*a*dm_dq;

% now, form jacobian

    jac_r = zeros(3,3);

    if ( abs(m) < 1 )
      jac_r(1,:) = dfa_dq;
      jac_r(2,:) = (dfa_dq*fb + dfb_dq*fa)*(1/g);
      jac_r(3,:) = (1/2)*(1/(g^2-1))*(dfa_dq*fb*fb + 2*fa*fb*dfb_dq);
    elseif ( m <= -1 )
      jac_r(1,:) = drho_dq*u   + rho*du_dq;
      jac_r(2,:) = drho_dq*u*u + 2*rho*u*du_dq + dp_dq;
      jac_r(3,:) = ( drhoet_dq + dp_dq )*u + (rhoet + p)*du_dq;
    end
end

function [jac_L,jac_R] = jac_flux_roe(qL,qR)
g = 1.4;
A = @(u,et) [        0                  1               0   ;
              -(3-g)*u^2/2          (3-g)*u            (g-1);
              (g-1)*u^3-g*u*et   g*et-3*(g-1)/2*u^2     g*u];
dcdp = @(rho,u)    [ 1      0         0;
                     u      rho       0;
                    u^2/2  rho*u    1/(g-1)];
%
rhoL = qL(1);
uL = qL(2);
pL = qL(3);
etL = pL/(rhoL*(g-1)) + 1/2*uL^2;
htL = etL + pL/rhoL;
%
rhoR = qR(1);
uR = qR(2);
pR = qR(3);
etR = pR/(rhoR*(g-1)) + 1/2*uR^2;
htR = etR + pR/rhoR;
%
% Roe-averaged variables at interface between i & j
R = sqrt(abs(rhoR/rhoL));
RoeDensity = R*rhoL;
RoeVelocity = (R*uR + uL)/(R+1);
sq_vel = RoeVelocity^2;
RoeEnthalpy = (R*htR + htL)/(R+1);
RoeSoundSpeed2 = (g-1)*(RoeEnthalpy-0.5*sq_vel);
if RoeSoundSpeed2 <0
    error('RoeSoundSpeed2 negative!');
end
RoeSoundSpeed = sqrt(RoeSoundSpeed2);
% Eigen values and Eigen vectors
LambdaVec  = [RoeVelocity, RoeVelocity+RoeSoundSpeed, RoeVelocity-RoeSoundSpeed];
Delta = 0.1;
for ii = 1:3
    LambdaVec(ii) = max(abs(LambdaVec(ii)), Delta*max(LambdaVec));
end
Lambda = diag(LambdaVec);
%
[P,Pinv] = GetPPinvMatrix(RoeDensity, RoeVelocity, RoeSoundSpeed);
%
jac_L = (0.5*A(uL,etL) + 0.5*P*Lambda*Pinv) * dcdp(rhoL,uL);
jac_R = (0.5*A(uR,etR) - 0.5*P*Lambda*Pinv) * dcdp(rhoR,uR);
end

function [P,Pinv] = GetPPinvMatrix(rho,u,a)
g = 1.4;
P = [1      rho/(2*a)                       -rho/(2*a)
     u      rho/(2*a)*(u+a)                 -rho/(2*a)*(u-a)
     u^2/2  rho/(2*a)*(u^2/2+u*a+a^2/(g-1))  -rho/(2*a)*(u^2/2-u*a+a^2/(g-1))];

Pinv = [  (2*a^2 - g*u^2 + u^2)/(2*a^2),       (u*(g - 1))/a^2,     -(g - 1)/a^2;
         -(u*(2*a + u - g*u))/(2*a*rho), (a + u - g*u)/(a*rho),  (g - 1)/(a*rho);
         -(u*(2*a - u + g*u))/(2*a*rho), (a - u + g*u)/(a*rho), -(g - 1)/(a*rho)];
end

% function [jac_L,jac_R] = flux_jac_central(q_L,q_R)
% g = 1.4;
% Jac = @(rho,u,p)   [u       rho                                     0;
%                     u^2     2*rho*u                                 1;
%                     u^3/2   (p/(g-1) + rho*u^2/2 + p) + rho*u^2     u + u/(g-1)];
% 
% jac_L = 1/2*Jac(q_L(1),q_L(2),q_L(3));
% jac_R = 1/2*Jac(q_R(1),q_R(2),q_R(3));
% 
% end

function source_jac = jac_source_1D(dadx)
source_jac = zeros(3,3);
source_jac(2,3) = dadx;
end

function a = speed_of_sound(p, rho)
g = 1.4;
a = sqrt(g*p/rho);
end