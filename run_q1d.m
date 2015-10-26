%% Code for Quasi 1-D Convergent Divergent Nozzle
% This is the code for calculating flow inside a C-D Nozzle
%% Set input values
close all;clear; clc;
tic;
set(0,'defaultlinelinewidth',2,'defaultaxesfontsize',25,'defaultaxesfontname','Times');
flux_scheme = 'flux_roe';%vanLeer'; % 'flux_central','flux_roe'
% Time stepping settings
N = 256;
step_size = 0; % percent value for FD
CSA_flag = 1;
plot_iterations = 1;
plot_solution = 1;
iterations = 2000;
isConverged = 0;
toler = 1e-10;
iterout = 500;
resinit = [1,1,1];
explicit_flag = 0;
%Grid parameters
set(0,'defaultlinelinewidth',2,'defaultaxesfontsize',17);
neq = 3;
% Joe's example
nozzle_length = 2;
area_throat = 0.2;
des_var_b = 0.4;
%
des_var_b = des_var_b + step_size/100*des_var_b;
area_func = @(x) area_throat+des_var_b*( 1 + sin(pi*(x-0.5)) );
dadx_func = @(x) des_var_b*pi*cos(pi*(x-0.5));
dadx_p_func = @(x) pi*cos(pi*(x-0.5));
area_p_func = @(x) ( 1 + sin(pi*(x-0.5)) );
% Anderson
% nozzle_length = 3;
%  area_throat = 1;
% area_func = @(x) area_throat+2.2*( x.^2 );
% dadx_func = @(x) 2.2*2*x;
% x_f = linspace(-nozzle_length/2,-0.1,N+1); % Only converging nozzle
x_f = linspace(-nozzle_length/2,nozzle_length/2,N+1); % position of faces
area_f = area_func(x_f);
x_cc = ( x_f(1:end-1) + x_f(2:end))/2; %x_cc(2:end+1)=x_cc; x_cc(end+1)=x_cc(end);
dx = x_f(2:end) - x_f(1:end-1); %dx(2:end+1)=dx; dx(end+1)=dx(end);
area_cc = area_func(x_cc); % area_cc(2:end+1)=area_cc; area_cc(end+1)=area_cc(end);
dadx_cc = dadx_func(x_cc); % dadx_cc(2:end+1)=dadx_cc; dadx_cc(end+1)=dadx_cc(end);
% dadx_p_cc = dadx_p_func(x_cc);
cell_vol = area_cc.*dx;
% [x_f',area_f']
% [area_cc',dadx_cc']
% Reference values
t0      = 600.0; % K
p0      = 300000.0; % Pa
pback   = -1.0; %0.998*p0; %-1.0;
rho_set = 3.4619176;
u_set   = 38.5459;
p_set   = -1.0;
if pback<0
    mref_end = 0.1;
    mref_mid = 1;
    CFL_min = 1;
    CFL_max = 1000;
    CFL_max_it = 100;
else
    mref_end = 0.1;
    mref_mid = 0.6;
    CFL_min = 1;
    CFL_max = 2;
    CFL_max_it = 100;
end
% Gas properties
r = 287; % J/Kg/K
gamma = 1.4;
%% Exact Solution
prim_exact = q1d_exact(t0,p0,x_cc,area_cc,area_throat,gamma,r,pback,area_f(end));
% load('prim_exact_data.mat');
% prim_exact_interp = interp1(prim_exact_data(:,1),prim_exact_data(:,2:4),x_cc,'pchip');
% prim_exact = prim_exact_interp.';
rhobar_exact = prim_exact(1,:)/(p0/(r*t0));
pbar_exact = prim_exact(3,:)/p0;
psi_exact = (1./pbar_exact).^((gamma-1)/(gamma));
M_exact = prim_exact(2,:)./sqrt(gamma*r*t0./psi_exact);
if plot_solution == 1
    fig_flow_soln=figure();
    plot(x_cc,area_cc','k');
    hold on;
    figure(fig_flow_soln);
    plot(x_cc,[rhobar_exact.',M_exact.',pbar_exact.']);
    % legend('Nozzle','\rho_{bar}','M','p_{bar}','location','nw');
    xlabel('x');title('Exact solution to Q-1D Nozzle');
    grid on; hold on;
end
%% Initialize Solution
% For isentropic expansion
prim_cc = zeros(neq,N);
% center = ceil(N/2);
m = zeros(1,length(x_cc));
for i = 1:length(x_cc)
    if pback < 0
        % Linear throughout  M variation for isentropic supersonic flow
        m(i) = mref_mid + x_cc(i)*(mref_mid - mref_end);
    else
        % Wedge M variation
        if x_cc(i)<=0
            m(i) = mref_mid + x_cc(i)*(mref_mid - mref_end);
        else
            m(i) = mref_mid + x_cc(i)*(mref_end - mref_mid);
        end
    end
    psi = 1 + 1/2*(gamma-1)*m(i)^2;
    t = t0/psi;
    p = p0/( psi^(gamma/(gamma-1)) );
    prim_cc(1,i) = p/(r*t);           % rho
    prim_cc(2,i) = m(i)*sqrt(gamma*r*t); % u
    prim_cc(3,i) = p;                 % p
end

% % Set initial condition as the exact solution-- JUST FOR CHECKING
% prim_cc = prim_exact;
%
scrsz = get(0,'ScreenSize');
% figure('Position',[scrsz(3)/6 scrsz(4)/3 4*scrsz(3)/5 scrsz(4)/2]);
% subplot(1,3,1)
% plot(x_cc,prim_cc(1,:)');title('Initial \rho');
% subplot(1,3,2)
% plot(x_cc,prim_cc(2,:)');title('Initial u');
% title('INITIAL CONDITIONS');
% subplot(1,3,3)
% plot(x_cc,prim_cc(3,:)');title('Initial p');
% % plot(x_cc,m);
% cons_cc = prim_to_cons(prim_cc); % assumed gamma = 1.4
%% Create output file headers
fp1 = fopen('./history.dat','w');
fprintf(fp1,'TITLE = "Q-1D Nozzle Iterative Residual History"\n');
fprintf(fp1,'variables="Iteration""Time(s)""Res1""Res2""Res3"\n');
fclose(fp1);
fp2 = fopen('./nozzle.dat','w');
fprintf(fp2,'TITLE = "Nozzle Field Data"\n');
fprintf(fp2,'variables="x""area""rho""u""p"\n');
fclose(fp2);
% Header for Screen Output
fprintf('Iter.   Time (s)       dtmin (s)         Continuity      x-Momentum     Energy\n');
%
write_output(0, x_cc, area_cc, prim_cc);
L2Mat = ones(3,iterations);
rtime = 0;
%
%% Set gc values explicitly only for the initial time step
prim_gi = subsonic_inflow_explicit(prim_cc(:,1), t0, p0);
prim_go = outflow_explicit(prim_cc(:,N), pback);
prim_gc = [prim_gi,prim_go];
%
for n = 1:iterations
    %% Set time step
    CFL = CFL_min + (CFL_max - CFL_min)/(CFL_max_it - 1)*(n-1);
    CFL = min(CFL,CFL_max);
    dt = set_time_step(dx, prim_cc,CFL);
    rtime = rtime + min(dt);
    
    %% Form Residual
    [Res, Res_gi, Res_go] = create_residual(flux_scheme, prim_cc, prim_gc, area_f, dadx_cc, dx, t0, p0);
    
    %% Check convergence
    [L2, resinit, conv] = check_iterative_convergence(n, Res, resinit, rtime, min(dt));
    L2Mat(:,n) = L2;
    if(conv<toler)
        fp1 = fopen('./history.dat','a');
        fprintf(fp1, '%d %e %e %e %e\n',n, rtime, L2(1), L2(2), L2(3));
        fclose(fp1);
        isConverged = 1;
        n_conv = n;
        break;
    end
    % Output solution file every 'iterout' steps
    if( (mod(n,iterout)==0) )
        write_output(n, x_cc, area_cc, prim_cc);
    end
    %% Form LHS & RHS
    if explicit_flag ~= 1
        [L, D, U, DGi, UGi, LGo, DGo] = fill_lhs(flux_scheme, dx, area_f, dadx_cc, prim_cc, prim_gc, pback);
        Jacobian = zeros(3*N+6,3*N+6);
        RHS = zeros(3*N+6,1);
        for row = 2:N+1
            cell = row-1;
            Jacobian(3*(row-1)+1:3*(row-1)+3 , 3*(row-2)+1:3*(row-2)+9) = ...
                                    [L(:,:,cell),D(:,:,cell),U(:,:,cell)];
            RHS(3*(row-1)+1:3*(row-1)+3,1) = -Res(:,cell);
        end
        % Add BC row blocks to the Jacobian and RHS
        Jacobian(1:3 , 1:6) = [DGi, UGi];
        RHS(1:3, 1) = -Res_gi;
        Jacobian(3*N+4:3*N+6 , 3*N+1:3*N+6) = [LGo, DGo];
        RHS(3*N+4:3*N+6, 1) = -Res_go;
        % Add time term
        LHS_Matrix = Jacobian;
        for row = 2:N+1
            cell = row-1;
            dcdp = [1                    0                               0;
                    prim_cc(2,cell)      prim_cc(1,cell)                 0;
                    prim_cc(2,cell)^2/2  prim_cc(1,cell)*prim_cc(2,cell) 1/(gamma-1)];
            %
            LHS_Matrix(3*(row-1)+1:3*(row-1)+3 , 3*(row-1)+1:3*(row-1)+3 ) = ...
              LHS_Matrix(3*(row-1)+1:3*(row-1)+3 , 3*(row-1)+1:3*(row-1)+3 ) + ...
              dcdp*cell_vol(cell)/dt(cell);
        end
    end
    
    %% Solve system of equations
    if explicit_flag == 1
        delta_prim_cc = zeros(3,N);
        for col = 1:N
            delta_prim_cc(:,col) = -Res(:,col)*dt(col);
        end
        % Update primitive variables and calculate ghost cell values
        prim_cc = prim_cc + delta_prim_cc;
        %
        prim_gi = subsonic_inflow_explicit(prim_cc(:,1), prim_cc(:,2), t0, p0);
        prim_go = outflow_explicit(prim_cc(:,N), prim_cc(:,N-1), pback);
        prim_gc = [prim_gi,prim_go];
    else
        RHS = RHS(:);
        delta_prim_cc_soln = LHS_Matrix\RHS;
        delta_prim_cc = reshape(delta_prim_cc_soln,3,N+2);
        % Update primitive variables and set ghost cell values
        prim_cc(:,1:N) = prim_cc(:,1:N) + delta_prim_cc(:,2:N+1);
        prim_gc = prim_gc + delta_prim_cc(:,[1,N+2]);
    end
    % Plot current solution
    if plot_iterations == 1
        if n==1
            scrsz = get(0,'ScreenSize');
            fig_current_sol = figure('Position',[scrsz(3)/6 scrsz(4)/3 4*scrsz(3)/5 scrsz(4)/2]);
            %             hold on;
        else
            rhobar = prim_cc(1,:)/(p0/(r*t0));
            pbar = prim_cc(3,:)/p0;
            psi = (1./pbar).^((gamma-1)/(gamma));
            M = prim_cc(2,:)./sqrt(gamma*r*t0./psi);
            figure(fig_current_sol)
            subplot(1,3,1)
            plot(x_cc,rhobar);
            legend(['Iter=',num2str(n)]);
            xlabel('x');ylabel('\rho');title('Density');
            grid on;
            subplot(1,3,2)
            plot(x_cc,M);
            xlabel('x');ylabel('M');title('Mach number');
            grid on;
            subplot(1,3,3)
            plot(x_cc,pbar);%title('p');
            xlabel('x');ylabel('p');title('Pressure');
            grid on;
        end
    end
        %
    
    %% Floor variables for stability
    for j = 1:N
        if prim_cc(1,j)<=0 || prim_cc(3,j)<=0
            prim_cc(1,j) = 0.01;
            prim_cc(2,j) = 1.0;
            prim_cc(3,j) = 3000.0;
        end
%         if prim_cc(2,j)>=100
%             prim_cc(2,j) = 100;
%         end

        if prim_cc(2,j)<=0
            prim_cc(2,j) = 50;
        end
    end
    
    %% Write solution
%     write_solution(n, N, x_cc, prim_cc, cons_cc);
%     
end
if isConverged ~= 1
    fprintf('Solution failed to converge in %d iterations!!!\n', iterations);
    n_conv = iterations;
end

if isConverged == 1
    fprintf('Solution converged in %d iterations!!!\n', n);
end
% Output solution and restart file
write_output(n, x_cc, area_cc, prim_cc);
fclose all;
%
if plot_solution == 1
    figure()
    semilogy(L2Mat(:,1:n_conv).')
    xlabel('Iteration number');ylabel('L_2 Norms');
    title('Flow solution congergence');
    legend('Cont.','Mom.','En.');
end
%
% figure()
% semilogy(Linf')
% xlabel('Iteration number');ylabel('L_{\infty} Norms');
% title('Congergence - L_{infty} norm');
% legend('\rho','u','p');
%
rhobar = prim_cc(1,:)/(p0/(r*t0));
pbar = prim_cc(3,:)/p0;
psi = (1./pbar).^((gamma-1)/(gamma));
M = prim_cc(2,:)./sqrt(gamma*r*t0./psi);
if plot_solution == 1
    figure(fig_flow_soln)
    plot(x_cc,[rhobar.',M.',pbar.'],'o');
    legend({'Nozzle','$\hat{\rho}$','$\hat{M}$','$\hat{p}$',...
        '$\rho$','$M$','$p$'},...
        'location','nw','Interpreter','latex');
    xlabel('x');title('CFD solution to Q-1D Nozzle');
    grid on;
end
% Plot T0 and P0 calculated values
M_cal = prim_cc(2,:)./sqrt(gamma*prim_cc(3,:)./prim_cc(1,:));
psi_cal = 1 + (gamma-1)/2*M_cal.^2;
T_cal = prim_cc(3,:)./(r*prim_cc(1,:));
T0_cal = T_cal.*psi_cal;
P0_cal = prim_cc(3,:).*psi.^(gamma/(gamma-1));
% figure()
% subplot(1,2,1)
% plot(x_cc, T0_cal);
% xlabel('x');title('T0 calculated');
% subplot(1,2,2)
% plot(x_cc, P0_cal);
% xlabel('x');title('P0 calculated');
%% Run CSA
save('Q1D_flow_soln.mat');
time_primary = toc;
fprintf('Primary analysis completed in %3.2f seconds!!!\n', time_primary);
if CSA_flag == 1
    run_q1d_CSA();
end