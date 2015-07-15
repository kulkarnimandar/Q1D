function prim_exact = q1d_exact(t0,p0,x_cc,area_cc,area_throat,g,R,p_back,area_back)

rho_0 = p0/(R*t0);
psi = @(M) 1 + (g-1)/2*M^2;
rho = @(M) rho_0/psi(M)^(1/(g-1));
u = @(M) M*sqrt(g*R*t0/psi(M));
p = @(M) p0/psi(M)^(g/(g-1));
%
prim_exact = zeros(3,length(x_cc));

if p_back<0
    Astar = area_throat;
else
    p_bar_inv = p0/p_back;
    psi_back = p_bar_inv^((g-1)/g);
    M_back = sqrt( (2/(g-1)*(psi_back-1) ) );
    A_ratio_back = sqrt( (1/M_back^2) * ( ( (2/(g+1))*( 1+(g-1)/2*M_back^2 ))^((g+1)/(g-1)) ) );
    Astar = area_back/A_ratio_back;
end

for i = 1:length(x_cc)
    %
    f = @(M) ( 2./(g+1).*( 1 + (g-1)./2.*M.^2) ).^((g+1)/(g-1)) - (area_cc(i)./Astar).^2.*M.^2;
    %     figure();ezplot(f,[0,3.5]);
    if p_back<0
        if x_cc(i)<=0
            M = fzero(f,[0.000001,1-.000001]);
        elseif x_cc(i)>=0.05
            M = fzero(f,3.05);
        else % x_cc(i)>=0.2
            M = fzero(f,[1+.000001,3]);
        end
    else
        if x_cc(i)<=0
            M = fzero(f,0.5);
        elseif x_cc(i)>=0.05
            M = fzero(f,0.5);
        else % x_cc(i)>=0.2
            M = fzero(f,0.5);
        end
    end
    prim_exact(:,i) = [rho(M);u(M);p(M)];
end


end