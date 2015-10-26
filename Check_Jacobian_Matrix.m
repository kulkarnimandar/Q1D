%% Check Jacobian using complex step
%
close all;clear all;clc;
load('Q1D_flow_soln.mat');
Jacobian_analytic = Jacobian(4:3*(N+1), 4:3*(N+1));
%
Jacobian_CS = zeros(3*N, 3*N);
q_0 = prim_cc;
[Res_0, Res_gi_0, Res_go_0] = create_residual(flux_scheme, q_0, prim_gc, area_f, dadx_cc, dx, t0, p0);
Res_0 = Res_0(:);
delta = 1e-16;
for i = 1:N*3
    switch_vector = zeros(3*N,1);
    switch_vector(i,1) = 1i;
    q_delta = q_0 + reshape(delta*switch_vector,3,N);
    [Res_delta, Res_gi_delta, Res_go_delta] = create_residual(flux_scheme, q_delta, prim_gc, area_f, dadx_cc, dx, t0, p0);
    Res_delta = Res_delta(:);
    Jacobian_CS(:,i) = imag(Res_delta)/delta;
end
%
Jacobian_diff = zeros(size(Jacobian_CS));
for i = 1:3*N
    for j = 1:3*N
        if Jacobian_analytic(i,j) ~= 0
            Jacobian_diff(i,j) = (Jacobian_CS(i,j) - Jacobian_analytic(i,j))*100/Jacobian_analytic(i,j);
        else
            Jacobian_diff(i,j) = (Jacobian_CS(i,j) - Jacobian_analytic(i,j))*100;
        end
    end
end
%
surf(Jacobian_diff)
shading interp
colorbar
caxis([-200 200])
 view([-80, 34])