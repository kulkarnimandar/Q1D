clc; clear; close all;

%% Calculate true solutions
NMat0 = [32,64,128,256,512,1024,2048];
% NMat0 = [32,64,128];
flux_scheme = 'flux_vanLeer';%roe';
dir_name = ['Grid_convergence_',flux_scheme];
if exist(dir_name,'dir') == 0
    mkdir(dir_name);
% else
%     error('Directory already exists!');
end
h = NMat0(end)./NMat0;
DE_2Norm_rho = zeros(size(h));
DE_2Norm_u = zeros(size(h));
DE_2Norm_p = zeros(size(h));
DE_InfNorm_rho = zeros(size(h));
DE_InfNorm_u = zeros(size(h));
DE_InfNorm_p = zeros(size(h));
for ii=1:length(NMat0)
    NMat = NMat0;    
    N = NMat(ii);
    if exist([dir_name,'\','Q1D_sens_local_soln_',num2str(N),'.mat'],'file') == 0
        close all; clc;
        run_q1d;
        movefile('Q1D_sens_local_soln.mat',[dir_name,'\','Q1D_sens_local_soln_',num2str(N),'.mat'])
        movefile('Q1D_flow_soln.mat',[dir_name,'\','Q1D_flow_soln_',num2str(N),'.mat'])
    else
        load([dir_name,'\','Q1D_sens_local_soln_',num2str(N),'.mat']);
    end
    h = NMat0(end)./NMat0;
    DE_2Norm_rho(ii) = sqrt( sum( abs( prim_cc_local(1,:) - prim_local_exact(1,:) ).^2  )/N );
    DE_2Norm_u(ii)   = sqrt( sum( abs( prim_cc_local(2,:) - prim_local_exact(2,:) ).^2  )/N );
    DE_2Norm_p(ii)   = sqrt( sum( abs( prim_cc_local(3,:) - prim_local_exact(3,:) ).^2  )/N ); 
    DE_InfNorm_rho(ii) = max( abs( prim_cc_local(1,:) - prim_local_exact(1,:) ) );
    DE_InfNorm_u(ii)   = max( abs( prim_cc_local(2,:) - prim_local_exact(2,:) ) );
    DE_InfNorm_p(ii)   = max( abs( prim_cc_local(3,:) - prim_local_exact(3,:) ) );
    save([dir_name,'\','Grid_conv_study_Q1D.mat']);
end

RateOfConv = figure('Position',[50 50 scrsz(3)/1.1 scrsz(4)/1.7]);
set(gcf,'defaultlinelinewidth',2,'defaultaxesfontsize',13)
subplot(1,2,1)
loglog(h,[DE_2Norm_rho',DE_2Norm_u',DE_2Norm_p'],'o-');
hold on;
loglog(h,[DE_InfNorm_rho',DE_InfNorm_u',DE_InfNorm_p'],'o--');
xlabel('Mesh refinement parameter, h');
ylabel('L^2 and L^{\infty} norms of error in \rho^{\prime}, u^{\prime}, p^{\prime}');
title('Local Derivative Discretization Errors Norms: \epsilon_u = u^{\prime}_h - u^{\prime}_{exact}');
legend('||\epsilon_{\rho}||_{2}','||\epsilon_u||_{2}','||\epsilon_p||_{2}',...
       '||\epsilon_{\rho}||_{\infty}','||\epsilon_u||_{\infty}','||\epsilon_p||_{\infty}',....
       'location','southeast');
set(gca,'XTick',[1,5,10,20,40,80,160])
xlim([0 300]);
set(gca,'XMinorTick','on','YMinorTick','on');
p_hat_2Norm_p = zeros(length(h)-1,1);
p_hat_InfNorm_p = zeros(length(h)-1,1);
p_hat_2Norm_u = zeros(length(h)-1,1);
p_hat_InfNorm_u = zeros(length(h)-1,1);
p_hat_2Norm_rho = zeros(length(h)-1,1);
p_hat_InfNorm_rho = zeros(length(h)-1,1);
for jj = 1:length(h)-1
    p_hat_2Norm_p(jj) = log(DE_2Norm_p(jj+1)/DE_2Norm_p(jj))/log(h(jj+1)/h(jj));
    p_hat_InfNorm_p(jj) = log(DE_InfNorm_p(jj+1)/DE_InfNorm_p(jj))/log(h(jj+1)/h(jj));
    p_hat_2Norm_u(jj) = log(DE_2Norm_u(jj+1)/DE_2Norm_u(jj))/log(h(jj+1)/h(jj));
    p_hat_InfNorm_u(jj) = log(DE_InfNorm_u(jj+1)/DE_InfNorm_u(jj))/log(h(jj+1)/h(jj));
    p_hat_2Norm_rho(jj) = log(DE_2Norm_rho(jj+1)/DE_2Norm_rho(jj))/log(h(jj+1)/h(jj));
    p_hat_InfNorm_rho(jj) = log(DE_InfNorm_rho(jj+1)/DE_InfNorm_rho(jj))/log(h(jj+1)/h(jj));
end
%
figure(RateOfConv)
% set(gcf,'defaultlinelinewidth',2,'defaultaxesfontsize',13)
subplot(1,2,2)
semilogx(h(2:end),[p_hat_2Norm_rho,p_hat_2Norm_u,p_hat_2Norm_p],'o-');
hold on;
semilogx(h(2:end),[p_hat_InfNorm_rho,p_hat_InfNorm_u,p_hat_InfNorm_p],'o--');
xlabel('Mesh refinement parameter, h');
ylabel('Rate of convergence');
title('Rate of convergence of \rho^{\prime}, u^{\prime}, p^{\prime}');
set(gca,'XTick',[1,5,10,20,40,80])
set(gca,'XMinorTick','on','YMinorTick','on');
axis([0 150 0 2]);
legend('||\epsilon_{\rho}||_{2}','||\epsilon_u||_{2}','||\epsilon_p||_{2}',...
       '||\epsilon_{\rho}||_{\infty}','||\epsilon_u||_{\infty}','||\epsilon_p||_{\infty}',....
       'location','southeast');