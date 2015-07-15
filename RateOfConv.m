close all;clear all;
load('Data_Multiple_Layers.mat');
%
scrsz = get(0,'ScreenSize');
set(0,'defaultlinelinewidth',2,'defaultaxesfontsize',25,'defaultaxesfontname','Times');
% RateOfConv = figure('Position',[50 50 scrsz(3)/1.1 scrsz(4)/2.2]);
RateOfConv = figure('Position',[50 50 scrsz(3)/1.1 scrsz(4)/1.7]);
set(gcf,'defaultlinelinewidth',2,'defaultaxesfontsize',13)
subplot(1,2,1)
% plot(h,[DE_InfNorm_u_analytic',DE_InfNorm_u_9Layer'],'o--');
plot(h,[DE_InfNorm_u_analytic',DE_InfNorm_u_1Layer',...
        DE_InfNorm_u_2Layer',DE_InfNorm_u_2Layer_complete2TS',...
        DE_InfNorm_u_6Layer',DE_InfNorm_u_9Layer'],'o--');
xlabel('Mesh refinement parameter, h');
ylabel('Infinity norms of error in p^{\prime},  | \epsilon |_\infty');
title('Error norms of local derivatives');
% legend('Analytic Gradients','9-Layer SGR','location','northwest');
legend('Analytic Gradients','1L, 1+ O','2L, 2+ O',...
       '2L, 2 O','6L 5 0','9L, 5 O','location','se');
xlim([1 6]);
set(gca,'XTick',[1,2,4,6,8,10],'xscale','log','yscale','log')
grid on;
%
figure(RateOfConv)
subplot(1,2,2)
% plot(h(1:end-1),[p_hat_InfNorm_u_analytic,p_hat_InfNorm_u_9Layer],'o--');
plot(h(1:end-1),[p_hat_InfNorm_u_analytic,p_hat_InfNorm_u_1Layer,...
                 p_hat_InfNorm_u_2Layer,p_hat_InfNorm_u_2Layer_complete2TS, ...
                 p_hat_InfNorm_u_6Layer,p_hat_InfNorm_u_9Layer],'o--');
xlabel('Mesh refinement parameter, h');
ylabel('Rate of convergence');
title('Rate of convergence/ Order of accuracy');
axis([1 6 0 3]);
legend('Analytic Gradients','1L, 1+ O','2L, 2+ O',...
       '2L, 2 O','6L 5 0','9L, 5 O','location','se');
set(gca,'XTick',[1,2,4,6,8,10],'xscale','log')
grid on;
%
%% Figure for Aviation 2015
figure()
plot(h,[DE_InfNorm_u_analytic',DE_InfNorm_u_1Layer',...
        DE_InfNorm_u_2Layer_complete2TS',DE_InfNorm_u_6Layer'],'o--');
xlabel('Mesh refinement parameter, h');
ylabel('| p^{\prime}_{CSA} - p^{\prime}_{MS} |_\infty');
title('Error norm of  \partialp/\partialL');
legend('Ana. BC','1L','2L','6L','location','se');
xlim([1 6]);
set(gca,'XTick',[1,2,4,6,8,10],'xscale','log','yscale','log')
grid on;
%
figure()
plot(h(1:end-1),[p_hat_InfNorm_u_analytic,p_hat_InfNorm_u_1Layer,...
                 p_hat_InfNorm_u_2Layer_complete2TS,p_hat_InfNorm_u_6Layer],'o--');
xlabel('Mesh refinement parameter, h');
ylabel('Rate of convergence');
title('Rate of convergence of \partialp/\partialL');
axis([1 6 0 4]);
legend('Ana. BC','1L SGR','2L SGR','6L SGR','location','se');
set(gca,'XTick',[1,2,4,6,8,10],'xscale','log')
grid on;