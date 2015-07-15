function [L2, resinit, conv] = check_iterative_convergence(n, res, resinit, rtime, dtmin, CSA_flag)

if nargin == 5
    CSA_flag = 0;
end

fsmall = 0.01;

Norm2 = [norm(res(1,:),2);
           norm(res(2,:),2);
           norm(res(3,:),2)]/size(res,2);
%      Linf(:,n) = [norm(Res(1,:),inf); norm(Res(2,:),inf); norm(Res(3,:),inf)]./Linf_ini;
%
L2 = zeros(3,1);
for k = 1:3
    if n <= 5
        resinit(k) = max(Norm2(k),resinit(k));
    end
    L2(k) = Norm2(k)/(resinit(k)+fsmall);
end

conv = max(L2);

% Write iterative residuals every 10 iterations
if ( (mod(n,1)==0)||(n==1) )
    if CSA_flag ~= 1
        fp1 = fopen('./history.dat','a');
    else
        fp1 = fopen('./history_local.dat','a');
    end
    fprintf(fp1, '%d %e %e %e %e\n',n, rtime, L2(1), L2(2), L2(3) );
    fclose(fp1);
    fprintf('%d   %e   %e   %e   %e   %e\n',n, rtime, dtmin, L2(1), L2(2), L2(3) );
end

% Write header for iterative residuals every 200 iterations
if ( (mod(n,200)==0)||(n==1) )
    fprintf('Iter.   Time (s)       dt (s)         Continuity      x-Momentum     Energy\n');
end

end