function cons_cc = prim_to_cons(prim_cc)


% prim_cc = [rho, u, p]';
% cons_cc = [rho, rho*u, rho*et]'
% rho*et = p/(gamma-1) + 1/2*rho*u^2;
gamma = 1.4;

rho = prim_cc(1,:);
u = prim_cc(2,:);
p = prim_cc(3,:);

cons_cc(1,:) = rho;              
cons_cc(2,:) = rho.*u;
cons_cc(3,:) = p/(gamma-1) + 1/2*rho.*u.^2;
end