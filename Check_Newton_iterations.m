%% Checking Newton type iterations for linear equation
% f = c*x0 + d;
close all;clear;clc;
c = 2;
d = -4; %x0 = -d/c = 2;
f = @(x) c*x.^1 + d;
figure(1); 
ezplot(f,[0,5]); hold on;
grid on;xlabel('x');ylabel('f');title('Zero of f'); 
% axis([0 5 -10 60]);
% plot(-d/c,f(-d/c),'rsq');
figure(1)
plot(fzero(f,3),f(fzero(f,3)),'rsq');
%
dt = 10;
res = 1;
tol = 1e-4;
% x = -5+randn;
x=100;
iter = 0;
%
% figure(1);
% plot((x - f(x)/c),f(x - f(x)/c),'c*'); hold on;
%
while abs(res)>tol && iter<100
    iter = iter + 1;
    res = f(x);
    %
    figure(1)
    plot(x,res,'go');hold on;
    chk_str = [iter,x,res]
    %
    xnew = x - res/(c + 1/dt);
    x = xnew;
end
