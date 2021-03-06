%% Quasi 1-D CSA: run_q1d_CSA()
CSA_flag = 1;
iterations_CSA = 50;
dt_CSA_scalar = 1;
conv = 1;
toler = 1e-10;
isConverged = 0;
dadx_p_cc = dadx_p_func(x_cc);
area_p_f = area_p_func(x_f);
area_p_cc = area_p_func(x_cc);
dt_primary = zeros(size(dt));
% dt_primary = dt;
%% Initialize sensitivity solution
% For isentropic expansion
prim_cc_local = 0*ones(neq,N);
% prim_cc_local = [prim_cc_local_delta1pc_central,prim_cc_local_delta1pc_central(:,end)];
% load('./Finite_Difference/Subsonic_Supersonic_Case/prim_cc_local_delta1pc_central.mat');
% prim_cc_local = prim_cc_local_delta1pc_central;
% for ii = 1:length(x_cc)
%     prim_cc_local(:,ii) = -1*[0.048278639525954;
%                               -1.172102369895277e+002;
%                               1.160738552859373e+004]*x_cc(ii);
%                            
% end
scrsz = get(0,'ScreenSize');
% figure('Position',[scrsz(3)/6 scrsz(4)/3 4*scrsz(3)/5 scrsz(4)/2]);
% subplot(1,3,1)
% plot(x_cc,prim_cc_local(1,:)');title('Initial \rho^{\prime}');
% subplot(1,3,2)
% plot(x_cc,prim_cc_local(2,:)');title('Initial u^{\prime}');
% title('INITIAL CONDITIONS -- LOCAL SENSITIVITY');
% subplot(1,3,3)
% plot(x_cc,prim_cc_local(3,:)');title('Initial p^{\prime}');
%% Create output file headers
fp1 = fopen('./history_local.dat','w');
fprintf(fp1,'TITLE = "Q-1D Nozzle Iterative Residual History LOCAL DERIVATIVE"\n');
fprintf(fp1,'variables="Iteration""Time(s)""Res1""Res2""Res3"\n');
fclose(fp1);
fp2 = fopen('./nozzle_local.dat','w');
fprintf(fp2,'TITLE = "Nozzle Field Local Sensitivity Data"\n');
fprintf(fp2,'variables="x""area""rho_prime""u_prime""p_prime"\n');
fclose(fp2);
% Header for Screen Output
fprintf('Iter.   Time (s)       dtmin (s)         Continuity      x-Momentum     Energy\n');
%
write_output(0, x_cc, area_cc, prim_cc_local, CSA_flag);
L2Mat = ones(3,iterations_CSA);
rtime = 0;
%
%% Set gc values explicitly only for the initial time step
prim_gi_local = subsonic_inflow_explicit(prim_cc_local(:,1), t0, p0, CSA_flag, prim_gc(:,1));
prim_go_local = outflow_explicit(prim_cc_local(:,N), pback, CSA_flag);
prim_gc_local = [prim_gi_local,prim_go_local];
%
LHS_Matrix_local = LHS_Matrix;
for n = 1:iterations_CSA
    %% Set time step
%     CFL = CFL_min + (CFL_max - CFL_min)/(CFL_max_it - 1)*(n-1);
%     CFL = min(CFL,CFL_max);
%     dt = set_time_step(dx, prim_cc,CFL);
    dt_CSA = dt_CSA_scalar*ones(size(area_cc));
    rtime = rtime + min(dt_CSA);
    
    %% Form Residual
    [Res, Res_gi, Res_go] = create_residual(flux_scheme, prim_cc_local, prim_gc_local, area_f, dadx_cc, dx, t0, p0, CSA_flag, prim_cc, prim_gc, dadx_p_cc, area_p_f);
    
    %% Check convergence
    [L2, resinit, conv] = check_iterative_convergence(n, Res, resinit, rtime, min(dt), CSA_flag);
    L2Mat(:,n) = L2;
    if(conv<toler)
        fp1 = fopen('./history_local.dat','a');
        fprintf(fp1, '%d %e %e %e %e\n',n, rtime, L2(1), L2(2), L2(3));
        fclose(fp1);
        isConverged = 1;
        n_conv = n;
        break;
    end
    % Output solution file every 'iterout' steps
    if( (mod(n,iterout)==0) )
        write_output(n, x_cc, area_cc, prim_cc_local, CSA_flag);
    end
    %% Form RHS and update LHS with new time term
    if explicit_flag ~= 1
        RHS = zeros(3*N+6,1);
        for row = 2:N+1
            cell = row-1;
            RHS(3*(row-1)+1:3*(row-1)+3,1) = -Res(:,cell);
            %
            if n == 1
                LHS_Matrix_local(3*(row-1)+1:3*(row-1)+3 , 3*(row-1)+1:3*(row-1)+3 ) = ...
                    LHS_Matrix(3*(row-1)+1:3*(row-1)+3 , 3*(row-1)+1:3*(row-1)+3 ) + ...
                    eye(3)/dt_CSA(cell)*dt_primary(cell);
            end
        end
        RHS(1:3, 1) = - Res_gi;
        RHS(3*N+4:3*N+6, 1) = - Res_go;
    end    
%% Solve system of equations
    delta_prim_cc_local_vector = LHS_Matrix_local\RHS;
%% Reshaping and updating local derivatives
    delta_prim_cc_local = reshape(delta_prim_cc_local_vector,3,N+2);
    % Update primitive variables and set ghost cell values
    prim_cc_local(:,1:N) = prim_cc_local(:,1:N) + delta_prim_cc_local(:,2:N+1);
    prim_gc_local = prim_gc_local + delta_prim_cc_local(:,[1,N+2]);
    %
    if plot_iterations == 1
        if n==1
            scrsz = get(0,'ScreenSize');
            fig_current_sens_sol = figure('Position',[scrsz(3)/6 scrsz(4)/3 4*scrsz(3)/5 scrsz(4)/2]);
            hold on;
            subplot(1,3,1)
            plot(x_cc,prim_cc_local(1,:)');title('\rho^{\prime}');
            subplot(1,3,2)
            plot(x_cc,prim_cc_local(2,:)');title('u^{\prime}');
            subplot(1,3,3)
            plot(x_cc,prim_cc_local(3,:)');title('p^{\prime}');
        else
            figure(fig_current_sens_sol)
            subplot(1,3,1)
            plot(x_cc,prim_cc_local(1,:)');title('\rho^{\prime}');
            legend(['Iter=',num2str(n)]);
            subplot(1,3,2)
            plot(x_cc,prim_cc_local(2,:)');title('u^{\prime}');
            subplot(1,3,3)
            plot(x_cc,prim_cc_local(3,:)');title('p^{\prime}');
        end
    end
    
%     %% Floor variables for stability
%     for j = 1:N
%         if prim_cc(1,j)<=0 || prim_cc(3,j)<=0
%             prim_cc(1,j) = 0.01;
%             prim_cc(2,j) = 1.0;
%             prim_cc(3,j) = 3000.0;
%         end
% %         if prim_cc(2,j)>=100
% %             prim_cc(2,j) = 100;
% %         end
% 
%         if prim_cc(2,j)<=0
%             prim_cc(2,j) = 0;
%         end
%     end
    
    %% Write solution
%     write_solution(n, N, x_cc, prim_cc, cons_cc);
%     
end
if isConverged ~= 1
    fprintf('CSA failed to converge in %d iterations!!!\n', iterations_CSA);
    n_conv = n;
end

if isConverged == 1
    fprintf('CSA converged in %d iterations!!!\n', n);
end
% Output solution and restart file
write_output(n, x_cc, area_cc, prim_cc_local, CSA_flag);
fclose all;
%
if plot_solution == 1
    figure()
    semilogy(L2Mat(:,1:n_conv).')
    xlabel('Iteration number');ylabel('L_2 Norms');
    title('CSA congergence');
    legend('Cont. CSE','Mom. CSE','En. CSE');
end
%
%% Plot local sensitivity FD
if plot_solution == 1
    if pback<0
%         load('./Finite_Difference/Subsonic_Supersonic_Case/flux_vanLeer_delta_1_0p1/prim_cc_local_delta1pc.mat');
        load('./Finite_Difference/Subsonic_Supersonic_Case/prim_cc_local_FD_central_best.mat');
        %     load('./Finite_Difference/Subsonic_Supersonic_Case/prim_cc_local_delta0p1pc_central.mat');
    else
        load('./Finite_Difference/Subsonic_Case/prim_cc_local_delta1pc_ConvNozzle.mat');
        load('./Finite_Difference/Subsonic_Case/prim_cc_local_delta1pc_central_ConvNozzle.mat');
        prim_cc_local_delta1pc = prim_cc_local_delta1pc_ConvNozzle;
        prim_cc_local_delta1pc_central = prim_cc_local_delta1pc_central_ConvNozzle;
    end
%     x_cc_saved_1 = linspace(x_cc(1),x_cc(end),size(prim_cc_local_delta1pc,2));
    x_cc_saved_2 = linspace(x_cc(1),x_cc(end),size(prim_cc_local_FD_central_best,2));
    %
    prim_local_exact = q1d_sensitivities_exact(t0,p0,area_cc,area_throat,gamma,r,area_p_cc,prim_exact);
    %
    fig_sens_soln = figure('units','normalized','outerposition',[0 0 1 1]);
    set(0,'defaultlinelinewidth',2,'defaultaxesfontsize',22,'defaultaxesfontname','Times');
    subplot(1,3,1)
%     plot(x_cc_saved_1,prim_cc_local_delta1pc(1,:)');
%     hold on;
    plot(x_cc_saved_2,prim_cc_local_FD_central_best(1,:)','r');
    hold on;
    plot(x_cc,prim_local_exact(1,:)','k');
    xlabel('x');ylabel('\partial\rho/\partialb');
    % plot(x_cc_saved,prim_cc_local_delta0p1pc_central(1,:)','k');
    %
    subplot(1,3,2)
%     plot(x_cc_saved_1,prim_cc_local_delta1pc(2,:)');
%     hold on;
    plot(x_cc_saved_2,prim_cc_local_FD_central_best(2,:)','r');
    hold on;
    plot(x_cc,prim_local_exact(2,:)','k');
    xlabel('x');ylabel('\partialu/\partialb');
    title('Local sensitivities');
    % plot(x_cc_saved,prim_cc_local_delta0p1pc_central(2,:)','k');
    %
    subplot(1,3,3)
%     plot(x_cc_saved_1,prim_cc_local_delta1pc(3,:)');
%     hold on;
    plot(x_cc_saved_2,prim_cc_local_FD_central_best(3,:)','r');
    hold on;
    plot(x_cc,prim_local_exact(3,:)','k');
    xlabel('x');ylabel('\partialp/\partialb');
    % plot(x_cc_saved,prim_cc_local_delta0p1pc_central(3,:)','k');
    %
    figure(fig_sens_soln)
    subplot(1,3,1)
    plot(x_cc,prim_cc_local(1,:)','go');
%     legend('Forward FD 1%','Central FD 1%','Exact','CSA','location','ne');
    legend('Central FD 1e-5%','Exact','CSA','location','ne');
    grid on;
    %
    subplot(1,3,2)
    plot(x_cc,prim_cc_local(2,:)','go');
    grid on;
    %
    subplot(1,3,3)
    plot(x_cc,prim_cc_local(3,:)','go');
    grid on;
    %
end
save('Q1D_sens_local_soln.mat');
time_CSA = toc -  time_primary;
fprintf('CSA completed in %3.2f seconds!!!\n', time_CSA);
fprintf('Time for primary + CSA = %3.2f seconds!!!\n', time_CSA + time_primary);