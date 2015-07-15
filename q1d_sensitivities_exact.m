function prim_local_exact = q1d_sensitivities_exact(T0,P0,area_cc,area_throat,g,R,area_p_cc,prim_exact)

% prim_exact = q1d_exact(T0,P0,x_cc,area_cc,area_throat,g,R,p_back,area_back);

rho = prim_exact(1,:);
u = prim_exact(2,:);
p = prim_exact(3,:);
%
T = p./(rho*R);
a = sqrt(g*p./rho);
M = u./a;
phi = 1 + (g-1)/2*M.^2;
%
n = (g+1)/(g-1);
C = (2/(g+1))^n*( -2./M.^3.*phi.^n  + ...
                   n*(g-1)./M.*phi.^(n-1));
M_prime = 2*area_cc.*area_p_cc./(area_throat^2*C);
%
phi_prime = (g-1)*M.*M_prime;
T_prime = -T0./phi.^2.*phi_prime;
p_prime = -p.^2./P0.*g./(g-1).*phi.^(g/(g-1) - 1).*phi_prime;
rho_prime = (1./(R.*T)).*( p_prime - R.*rho.*T_prime );
a_prime = (1./(2.*a)).*( g.*p_prime./rho - g.*p./rho.^2.*rho_prime );
u_prime = a.*( M_prime + u./a.^2.*a_prime );
%
prim_local_exact = [rho_prime; u_prime; p_prime];
end