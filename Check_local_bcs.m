g = 1.4;
T0      = 600.0; % K
P0      = 300000.0; % Pa
R = 287; % J/Kg/K

% Values near the left boundary, 1st cell
rho = prim_cc(1,1);
u = prim_cc(2,1);
p = prim_cc(3,1);
T = rho/(g*R);

rho_p = prim_cc_local_delta1pc_central_ConvNozzle(1,1);
u_p = prim_cc_local_delta1pc_central_ConvNozzle(2,1);
p_p = prim_cc_local_delta1pc_central_ConvNozzle(3,1);
T_p = rho_p/(g*R);

a = sqrt(g*p/rho);
M = u/a;
a_p = (1/(2*a))*( g*p_p/rho - g*p*rho_p/rho^2);
M_p = u_p/a - u*a_p/a^2;

T0_p_by_T =  T0*T_p/T^2 + (g-1)*M*M_p 
P0_p_by_p =  P0*p_p/p^2 + g*(P0/p)*(T/T0)*M_p
