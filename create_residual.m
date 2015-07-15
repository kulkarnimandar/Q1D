function [Res, Res_gi, Res_go] = create_residual(prim_cc, prim_gc, area_f, dadx_cc, dx, T_stag_BC, P_stag_BC, CSA_flag, prim_cc_flow, prim_gc_flow, dadx_p_cc, area_p_f)
%
if nargin == 7
    CSA_flag = 0;
    prim_cc_flow = zeros(size(prim_cc));
    prim_gc_flow = zeros(size(prim_gc));
    dadx_p_cc = zeros(size(dadx_cc));
    area_p_f = zeros(size(area_f));
end
%
N = length(dadx_cc);
faces = N+1;
Res = zeros(3,N);
prim_gi = prim_gc(:,1);
prim_go = prim_gc(:,2);
%%% To check code: get exact values of q at the interfaces
% load('prim_exact_data.mat');
% prim_exact_interp = interp1(prim_exact_data(:,1),prim_exact_data(:,2:4),x_f,'pchip');
% prim_f = prim_exact_interp.';
%%%
% figure(10);title('Flux');xlabel('node #');ylabel('Flux');hold on;
ctr = 0;
%% Fill up residuals from cell #1 to cell #N
% Source terms added only for cells #1 to #(N-1)
for f=2:faces-1
    ctr = ctr+1;
    c=f-1;
    Flux = flux_vanLeer(prim_cc(:,c),prim_cc(:,c+1),CSA_flag,prim_cc_flow(:,c),prim_cc_flow(:,c+1));
    Res(:,c) = Res(:,c) + Flux*area_f(f);
    Res(:,c+1) = Res(:,c+1) - Flux*area_f(f);
    if CSA_flag == 1
        Flux_primary = flux_vanLeer(prim_cc_flow(:,c),prim_cc_flow(:,c+1));
        Res(:,c) = Res(:,c) + Flux_primary*area_p_f(f);
        Res(:,c+1) = Res(:,c+1) - Flux_primary*area_p_f(f);
    end
    % Add source term
    S = source(prim_cc(:,c),dadx_cc(c),CSA_flag,prim_cc_flow(:,c),dadx_p_cc(c));
    Res(:,c) = Res(:,c) - S*dx(c);
end
% legend('Central','van Leer');
%% Contribution of f1 to Res(:,1)
Flux = flux_vanLeer(prim_gi,prim_cc(:,1),CSA_flag,prim_gc_flow(:,1),prim_cc_flow(:,1));
Res(:,1) = Res(:,1) - Flux*area_f(1);
if CSA_flag == 1
    Flux_primary = flux_vanLeer(prim_gc_flow(:,1),prim_cc_flow(:,1));
    Res(:,1) = Res(:,1) - Flux_primary*area_p_f(1);
end
%% Contribution of f_(N+1) to Res(:,N)
Flux = flux_vanLeer(prim_cc(:,N),prim_go,CSA_flag,prim_cc_flow(:,N),prim_gc_flow(:,2));
Res(:,N) = Res(:,N) + Flux*area_f(N+1);
if CSA_flag == 1
    Flux_primary = flux_vanLeer(prim_cc_flow(:,N),prim_gc_flow(:,2));
    Res(:,N) = Res(:,N) + Flux_primary*area_p_f(N+1);
end
%% Add source term for the cell #N
S = source(prim_cc(:,N),dadx_cc(N),CSA_flag,prim_cc_flow(:,N),dadx_p_cc(N));
Res(:,N) = Res(:,N) - S*dx(N);
%
%% BC function value at the input and output ghost cell
% {psi_1, psi_2, psi_3}^T 
g = 1.4;
R = 287;
%
if CSA_flag ~= 1
    rho0 = prim_gi(1);
    u0 = prim_gi(2);
    p0 = prim_gi(3);
    u1 = prim_cc(2,1);
    %
    T0 = p0/(rho0*R);
    a0 = speed_of_sound(p0, rho0);
    M0 = u0/a0;
    phi = 1 + (g-1)/2*M0^2;
    T_stag = T0*phi;
    psi_1 = T_stag - T_stag_BC;
    %
    psi_2 = u0 - u1;
    %
    P_stag = p0*phi^(g/(g-1));
    psi_3 = P_stag - P_stag_BC;
else    
    rho0 = prim_gc_flow(1,1);
    u0 = prim_gc_flow(2,1);
    p0 = prim_gc_flow(3,1);
    %
    T0 = p0/(rho0*R);
    a0 = speed_of_sound(p0, rho0);
    M0 = u0/a0;
    phi = 1 + (g-1)/2*M0^2;
    %
    rho0_prime = prim_gi(1);
    u0_prime = prim_gi(2);
    p0_prime = prim_gi(3);
    u1_prime = prim_cc(2,1);
    %
    a0_prime = 1/(2*a0)*(g*p0_prime/rho0 - g*p0/rho0^2*rho0_prime);
    M0_prime = u0_prime/a0 - u0/a0^2*a0_prime;
    phi_prime = (g-1)/2*2*M0*M0_prime;
    T0_prime = p0_prime/(rho0*R) - p0/(rho0^2*R)*rho0_prime;
    %
    % T_stag = T0*phi;
    T_stag_prime = T0_prime*phi + T0*phi_prime;
    psi_1 = T_stag_prime - 0;
    %
    psi_2 = u0_prime - u1_prime;
    %
    % P_stag = p0*phi^(g/(g-1));
    P_stag_prime = p0_prime*phi^(g/(g-1)) + p0*(g/(g-1))*phi^((g/(g-1)) - 1)*phi_prime;
    psi_3 = P_stag_prime - 0;
end
%
Res_gi = [psi_1; psi_2; psi_3];
%
Res_go = prim_go - prim_cc(:,N);
%% Plot residuals
% Flux_throat = flux_central(prim_f(:,ceil(size(Res,2)/2)),prim_f(:,ceil(size(Res,2)/2)));
% scrsz = get(0,'ScreenSize');
% figure('Position',[scrsz(3)/6 scrsz(4)/3 4*scrsz(3)/5 scrsz(4)/2]);
% subplot(1,3,1)
% plot(Res(1,:)'/f1(1));title('Res Cont.');
% subplot(1,3,2)
% plot(Res(2,:)'/f1(2));title('Res Mom.');
% subplot(1,3,3)
% plot(Res(3,:)'/f1(3));title('Res Eng.');
end

function [Flux] = flux_vanLeer(qL,qR,CSA_flag,qL_flow,qR_flow)
g = 1.4;

if nargin == 2
    CSA_flag = 0;
%     qL_flow = zeros(size(qL));
%     qR_flow = zeros(size(qR));
end

if CSA_flag ~= 1
    % Calculate Left (+) Flux
    a = speed_of_sound(qL(3), qL(1));
    M = qL(2)/a;
    FL = zeros(3,1);
    %
    if ( abs(M) < 1 ) % Left subsonic flux
        fa = (1/4)*qL(1)*a*(M+1)^2;
        fb = a*((g-1)*M + 2);
        
        FL(1) = fa;
        FL(2) = fa*fb*(1/g);
        FL(3) = (1/2)*fa*fb*fb*(1/(g^2-1));
    elseif ( M >= 1 ) % Left supersonic flux
        FL(1) = qL(1)*a*M;
        FL(2) = qL(1)*a^2*(M^2+(1/g));
        FL(3) = qL(1)*a^3*M*((1/2)*M^2+(1/(g-1)));
    else
        FL = 0;
    end
    %
    % Calculate Right (-) Fluxes
    a = speed_of_sound(qR(3), qR(1));
    M = qR(2)/a;
    FR = zeros(3,1);
    %
    if ( abs(M) < 1 ) % Right subsonic flux
        fa = -(1/4)*qR(1)*a*(M-1)^2;
        fb = a*((g-1)*M - 2);
        
        FR(1) = fa;
        FR(2) = fa*fb*(1/g);
        FR(3) = (1/2)*fa*fb*fb*(1/(g^2-1));
    elseif ( M <= -1 ) % Right supersonic flux
        FR(1) = qR(1)*a*M;
        FR(2) = qR(1)*a^2*(M^2+(1/g));
        FR(3) = qR(1)*a^3*M*((1/2)*M^2+(1/(g-1)));
    else
        FR = 0;
    end
    %
    Flux = FL + FR;    
    
else % Case for CSA -- flux vector is F'
    
    
    % Flux calculated by central scheme
    q_face = (qL_flow + qR_flow)/2;
    q_p_face = (qL + qR)/2;
%     q_face = qL_flow;
%     q_p_face = qL;
    %
    rho = q_face(1);
    u = q_face(2);
    p = q_face(3);
    %
    rho_p = q_p_face(1);
    u_p = q_p_face(2);
    p_p = q_p_face(3);
    %
    %
    % % Flux vector as function of primitive variables for primary analysis
    %     Flux =   [rho*u;
    %               p+rho*u^2;
    %               u*p*g/(g-1) + 1/2*rho*u^3];
    %
    Flux_central =   [rho_p*u + rho*u_p;
        p_p + rho_p*u^2 + rho*2*u*u_p;
        u_p*p*g/(g-1) + u*p_p*g/(g-1) + 1/2*rho_p*u^3 + 1/2*rho*3*u^2*u_p];

    % Flux calculated by van Leer scheme
    % Calculate Left (+) Flux
    rho = qL_flow(1);
    u = qL_flow(2);
    p = qL_flow(3);
    %
    rho_p = qL(1);
    u_p = qL(2);
    p_p = qL(3);
    %
    a = speed_of_sound(p, rho);
    M = u/a;
    a_p = (1/(2*a))*( g*p_p/rho - g*p*rho_p/rho^2);
    M_p = u_p/a - u*a_p/a^2;
    %
    FL = zeros(3,1);
    %
    if ( abs(M) < 1 ) % Left subsonic flux
        fa = (1/4)*rho*a*(M+1)^2;
        fb = a*((g-1)*M + 2);
        %
        fa_p = (1/4)*rho_p*a*(M+1)^2 + (1/4)*rho*a_p*(M+1)^2 + (1/4)*rho*a*2*(M+1)*M_p;
        fb_p = a_p*((g-1)*M + 2) + a*(g-1)*M_p;
        %
        FL(1) = fa_p;
        FL(2) = fa_p*fb*(1/g) + fa*fb_p*(1/g);
        FL(3) = (1/2)*fa_p*fb*fb*(1/(g^2-1)) + (1/2)*fa*2*fb*fb_p*(1/(g^2-1));
    elseif ( M >= 1 ) % Left supersonic flux
        FL(1) = rho_p*a*M + rho*a_p*M + rho*a*M_p;
        FL(2) = rho_p*a^2*(M^2+(1/g)) + rho*2*a*a_p*(M^2+(1/g)) + rho*a^2*2*M*M_p;
        FL(3) = rho_p*a^3*((1/2)*M^3 + (M/(g-1))) + rho*3*a^2*a_p*((1/2)*M^3 + (M/(g-1))) + ...
            rho*a^3*((1/2)*3*M^2*M_p + (M_p/(g-1)));
    else
        FL = 0;
    end
    %
    % Calculate Right (-) Fluxes
    rho = qR_flow(1);
    u = qR_flow(2);
    p = qR_flow(3);
    %
    rho_p = qR(1);
    u_p = qR(2);
    p_p = qR(3);
    %
    a = speed_of_sound(p, rho);
    M = u/a;
    a_p = (1/(2*a))*( g*p_p/rho - g*p*rho_p/rho^2);
    M_p = u_p/a - u*a_p/a^2;
    %
    FR = zeros(3,1);
    %
    if ( abs(M) < 1 ) % Right subsonic flux
        fa = -(1/4)*rho*a*(M-1)^2;
        fb = a*((g-1)*M - 2);
        %
        fa_p = -(1/4)*rho_p*a*(M-1)^2 -(1/4)*rho*a_p*(M-1)^2 -(1/4)*rho*a*2*(M-1)*M_p;
        fb_p = a_p*((g-1)*M - 2) + a*(g-1)*M_p;
        %
        FR(1) = fa_p;
        FR(2) = fa_p*fb*(1/g) + fa*fb_p*(1/g);
        FR(3) = (1/2)*fa_p*fb*fb*(1/(g^2-1)) + (1/2)*fa*2*fb*fb_p*(1/(g^2-1));
    elseif ( M <= -1 ) % Right supersonic flux
        FR(1) = rho_p*a*M + rho*a_p*M + rho*a*M_p;
        FR(2) = rho_p*a^2*(M^2+(1/g)) + rho*2*a*a_p*(M^2+(1/g)) + rho*a^2*2*M*M_p;
        FR(3) = rho_p*a^3*((1/2)*M^3 + (M/(g-1))) + rho*3*a^2*a_p*((1/2)*M^3 + (M/(g-1))) + ...
            rho*a^3*((1/2)*3*M^2*M_p + (M_p/(g-1)));
    else
        FR = 0;
    end
    Flux = FL + FR;
    Flux_diff = Flux - Flux_central;
    
end

end

function [Flux] = flux_central(q_L,q_R)
g = 1.4;
% rho_et = @(rho,u,p) p/(g-1) + 1/2*rho*u^2;
% F = @(rho,u,p)   [rho*u;
%                   p+rho*u^2;
%                   u*(rho_et(rho,u,p) + p)];
% 
% Flux = F( 0.5*(q_L(1)+q_R(1)) , 0.5*(q_L(2)+q_R(2)) , 0.5*(q_L(3)+q_R(3)) );

F_cons = @(Q1,Q2,Q3) [Q2;
                      Q2^2/Q1 + (g-1)*(Q3 - 1/2*Q2^2/Q1);
                      Q2*Q3/Q1 + (g-1)*Q2/Q1*(Q3 - 1/2*Q2^2/Q1)];
%                            
Q_L = prim_to_cons([q_L(1);q_L(2);q_L(3)]);
Q_R = prim_to_cons([q_R(1);q_R(2);q_R(3)]);
Flux = F_cons( 0.5*(Q_L(1)+Q_R(1)) , 0.5*(Q_L(2)+Q_R(2)) , 0.5*(Q_L(3)+Q_R(3)) );

end

function S = source(q,dadx,CSA_flag,q_flow,dadx_p)
if CSA_flag ~= 1
    S = zeros(3,1);
    S(2) = q(3)*dadx;
else
    S = zeros(3,1);
    p = q_flow(3);
    p_p = q(3);
    S(2) = p_p*dadx + p*dadx_p;
end
end

function a = speed_of_sound(p, rho)
g = 1.4;
a = sqrt(g*p/rho);
end