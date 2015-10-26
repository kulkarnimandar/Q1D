%% Q1D Finite Difference Step size study
close all;clear; clc;
CSA_flag = 0;
percMat = [0,logspace(-10,0,101)];
N = 256;
flux_scheme = 'flux_vanLeer';%roe';
plot_iterations = 0;
plot_solution = 0;
dir_name = ['Finite_Difference\Subsonic_Supersonic_case\',flux_scheme];
if exist(dir_name,'dir') == 0
    mkdir(dir_name);
else
    error('Directory already exists!');
end
Sens_2Norm_FD_fwd = zeros(3,length(percMat));
Sens_2Norm_FD_central = zeros(3,length(percMat));
for ii=1:length(percMat)
    step_size = percMat(ii);
    fprintf('Running Q1D for step size %3.6f percent\n',step_size);
    run_q1d;
    if ii == 1
        prim_cc_0 = prim_cc;
        des_var_b_0 = des_var_b;
        continue;
    else
        prim_cc_delta = prim_cc;
        delta_b = step_size/100*des_var_b_0;
        prim_cc_local_FD_fwd = (prim_cc_delta - prim_cc_0)/delta_b;
    end
    %
    step_size = -percMat(ii);
    fprintf('Running Q1D for step size %3.6f percent\n',step_size);
    run_q1d;
    prim_cc_minus_delta = prim_cc;
    prim_cc_local_FD_central = (prim_cc_delta - prim_cc_minus_delta)/(2*delta_b);
    %
    for jj = 1:3
        Sens_2Norm_FD_fwd(jj,ii) = norm(prim_cc_local_FD_fwd(jj,:),2);
        Sens_2Norm_FD_central(jj,ii) = norm(prim_cc_local_FD_central(jj,:),2);
    end
    save([dir_name,'\','Q1D_step_size_study.mat']);
end
%
set(0,'defaultlinelinewidth',2,'defaultaxesfontsize',13,'defaultaxesfontname','Times');
figure()
loglog(percMat,Sens_2Norm_FD_fwd(3,:),'r-');
hold on;
loglog(percMat,Sens_2Norm_FD_central(3,:),'k-');
xlabel('Step size \Delta b');%,'interpreter','latex');
ylabel('$||p^{\prime}||_2\,\,\,$','interpreter','latex');
set(gca, 'xdir','reverse')
title('FD step size study');
legend('Forward finite difference','Central finite difference',...
       'location','sw');
axis([1e-10 1e0 10^5.84749 10^5.84751]);
%
% Best step size = 1e-5
step_size_best = 1e-5;
step_size = step_size_best;
fprintf('Running Q1D for step size %3.6f percent\n',step_size);
run_q1d;
prim_cc_delta = prim_cc;
delta_b = step_size/100*des_var_b_0;
step_size = -step_size_best;
fprintf('Running Q1D for step size %3.6f percent\n',step_size);
run_q1d;
prim_cc_minus_delta = prim_cc;
prim_cc_local_FD_central_best = (prim_cc_delta - prim_cc_minus_delta)/(2*delta_b);
%
fig_sens_soln = figure('units','normalized','outerposition',[0 0 1 1]);
set(0,'defaultlinelinewidth',2,'defaultaxesfontsize',22,'defaultaxesfontname','Times');
subplot(1,3,1)
plot(x_cc,prim_cc_local_FD_central_best(1,:)','b');
xlabel('x');ylabel('\partial\rho/\partialb');
% plot(x_cc_saved,prim_cc_local_delta0p1pc_central(1,:)','k');
%
subplot(1,3,2)
plot(x_cc,prim_cc_local_FD_central_best(2,:)','b');
xlabel('x');ylabel('\partialu/\partialb');
title('Local sensitivities');
legend(['Best FD step size = ',num2str(step_size_best)]);
%
subplot(1,3,3)
plot(x_cc,prim_cc_local_FD_central_best(3,:)','b');
xlabel('x');ylabel('\partialp/\partialb');
save([dir_name,'\','Q1D_step_size_study.mat']);