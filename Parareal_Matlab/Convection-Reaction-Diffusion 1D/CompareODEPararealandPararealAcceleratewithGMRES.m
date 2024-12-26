close all
clear
clc


%% ---------------------------------------------------------------
%  Parareal algorithm for linear problem:
%         y'(t) = a*y(t), t in (0,T)
%          y(t=0) = y0


% non-overlapping subdomains 2-level domain decomposition  preconditioner

global a N m K T;
a = -1;
yExactFunction = @(t)exp(a*t);

% stability function (  y^n+1 = R(z)y^n  )
Rz = @(z) 1/(1-z);     % Backward Euler
% Rz = @(z) 1+z;       % Forward Euler

%% -----------------------------------------------------------------


% initial condition
u0 = 1;

% coarse propagator
N = 100;     %   number of coarse time intervals
T = 100;      %   T = N*Dt
Dt = T/N;   %   coarse timesteps

xt = linspace(0,T,N+1);

% fine propagator
m = 20;    % number of fine time steps in each coarse time interval
dt = Dt/m;

K = 20 % number of Parareal iterations

cmap = hsv(K); % color map

% 2-level domain decomposition initial
N_fine_nodes = m*N+1; % 11
% N_subdomains = m+1; % 6
Id = eye(N_fine_nodes);

% plot the time domain
figure
plot(0:N_fine_nodes-1,1,'.b',0:m:N_fine_nodes-1,1,'ro')
xlim([0 N_fine_nodes-1])
ylim([0.5 1.5])

%%
% original problem matrix
phi = (1-a*dt)^-1;
A = eye(m*N+1,m*N+1) + (-phi)*diag(ones(m*N,1),-1);

% Reduced system

% fine propagator
F_tilde = (1-a*dt)^-m;

% coarse propagator
G_tilde = (1-a*Dt)^-1;



% Uk_0 = zeros(m*N,1);
% Uk_0(1) = u0;

rhs0 = [u0; zeros(m*N,1)];
% rhs = [u0; phi*u0;-G_tilde*u0;phi*U(3);-G_tilde*U(3)];


% exact solution

ttrue = linspace(0,T,N*m+1);
pick = 1:m:N*m+1;

uF_true = A\rhs0;

% reduced matrix
A_tilde = eye(N+1,N+1) + (-F_tilde)*diag(ones(N,1),-1);
% A_tilde = eye(m*N+1,m*N+1) + (-F_tilde)*diag(ones(m*N,1),-1);

% approximate reduced matrix
B = eye(N+1,N+1) + (-G_tilde)*diag(ones(N,1),-1);




% coarse space initial
indice_R_coarse  = 0:m:m*N;
R = zeros(length(indice_R_coarse),N_fine_nodes);
for i = 1:length(indice_R_coarse)
    R(i,indice_R_coarse(i)+1) = 1;
end

% first subdomain  Omega_1 = {0}
R1 = zeros(1,N_fine_nodes);
R1(1) = 1;
A1 = R1*A*R1';

% define nodes in other subdomains
N_nodes_subdomain = m;
indice_R_subdomains = zeros(N,N_nodes_subdomain);
for i =1:N
    indice_R_subdomains(i,:) =  (i-1)*m+1:(i-1)*m+m;
end

% Omega_2 = {1,2} Omega_3 = {3,4}
R2_end = zeros(N*N_nodes_subdomain,N_fine_nodes) ;
for i = 1:N
    for j = 1:N_nodes_subdomain
        R2_end(N_nodes_subdomain*(i-1)+j,indice_R_subdomains(i)+j) = 1;
    end
end

R_sub = [];       % restriction matrices
A_sub = [];       % sub-matrices
for i = 2: N+1
    R_sub{i} = R2_end(N_nodes_subdomain*(i-2)+1:N_nodes_subdomain*(i-2)+N_nodes_subdomain,:);
    A_sub{i} = R_sub{i}*A*R_sub{i}';
end

% define 2-level domain decomposition multiplicative preconditioner
P = R1'*A1*R1;
for i = 2:N+1
    P = P + R_sub{i}'*A_sub{i}*R_sub{i};
end

%Multiplicative
M = P*(R'*B*R +  Id-R'*R);

%Additive
% M = P + (R'*B*R - R'*R)

% preconditioned stationary iteration
% U^(k+1) = U^(k) + M^(-1)*( f-A*U^(k)) )




%% ------------------------------------------------------------------

% initial coarse integration

B1 =  eye(N+1,N+1) + (-G_tilde)*diag(ones(N,1),-1);
% B1 =  eye(N*m+1,N*m+1) + (-G_tilde)*diag(ones(N*m,1),-1);

%Initial coarse Parareal iteration at the coarse level
rhs0_coarse = [u0;zeros(N,1)];
U0_coarse = B1\rhs0_coarse;

U_fine=spline(ttrue(pick),U0_coarse,ttrue)';
U0_fine = U_fine;
% first 2-level dd preconditioning iterartion

%Initial fine solution at the fine level

% rhs=[u0];
% for i=1:N
%     rhs=[rhs; phi*U0_coarse(i);-G_tilde*U0_coarse(i)];
% end
% U_fine = M\rhs;

% U_fine=zeros(N_fine_nodes,1);
% U_fine(1:2:end,1) =  U0_coarse;
% for i=1:2:N_fine_nodes-2
%     U_fine(i+1) = F_tilde*U_fine(i);
% end

% rhs = zeros(N_fine_nodes,1);
% rhs(1) = u0;
% U_fine = A_tilde\rhs;

figure
hold on
plot(ttrue,uF_true,'-^k',ttrue,U_fine,'-*m');
legend('True solution','Parareal solution at k = 0')

U0n = U_fine;
proNj = 1;

L2NormError = sqrt(sum(Dt*(uF_true(pick)-U_fine(pick)).^2));%norm(uF_true(pick)-U_fine(pick),2);%
LInfNormError = norm(uF_true(pick)-U_fine(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));
SupperlinearErrorBound  = LInfNormError;
LinearErrorBound  = LInfNormError;

%%
% rhs
f = zeros(N_fine_nodes,1);
f(1) = u0;


% norm(f-A*U0_fine)
% norm(f-A*U_fine)
res_Parareal = 1; % relative residual ||r||/||r0||
% res_Parareal(1,:) = norm(f-A*U0_fine);
% iteration loop
for k = 1: K
    % k = k + 1
% residual
    r_k = f - A*U_fine;
    
  % first level ~ fine propagator ~ additive Schwarz in parallel
%     U_fine_temp = P\r_k;   
    for i = 1:N+1
%         i = i + 1
        if i == 1
            x_temp_1 = R1*r_k;
            x_temp_2 = A1\x_temp_1;
            x{i} = R1'*x_temp_2;
            U_fine_temp = x{i};         
        else
            x_temp_1 = R_sub{i}*r_k;
            x_temp_2 = A_sub{i}\x_temp_1;
            x{i} = R_sub{i}'*x_temp_2;  
            U_fine_temp = U_fine_temp + x{i};
        end            
    end

  % second level ~ coarse grid correction sequentially
    U_fine_temp_1 = R*U_fine_temp;
    U_fine_temp_2 = B\U_fine_temp_1;
    U_fine_temp_3 = R'*U_fine_temp_2;
    U_fine_temp_4 = U_fine_temp_3 + U_fine_temp - R'*R*U_fine_temp;
    U_fine = U_fine + U_fine_temp_4;
    
%     U_fine = U_fine + (R'*M_tilde*R + Id - R'*R)\U_fine_temp;
    
%     U_fine = M\(M*U_fine + (f-A*U_fine));
    res_Parareal(k+1,:) = norm(f-A*U_fine)/norm(f-A*U0_fine);
    
%      U_fine = M\(M*U_fine + (f-A*U_fine));
    %         U_fine = gmres(M,M*U_fine + (f-A*U_fine),10,1e-6);
    
    plot(ttrue(pick),U_fine(pick),'Color',cmap(k,:),'Marker','o');
    
    %
    %        norm(uF_true(pick)-U,'inf')
    proNj(k+1) = proNj(k)*(N-k);

    L2NormError(k+1) = sqrt(sum(Dt*(uF_true(pick)-U_fine(pick)).^2));%norm(uF_true(pick)-U_fine(pick),2);%
    LInfNormError(k+1) = norm(uF_true(pick)-U_fine(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));%
    SupperlinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt)))^(k)/factorial(k)*proNj(k+1)*max(abs(uF_true(pick)-U0n(pick)));
    LinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt))/(1- abs(Rz(a*Dt))))^(k)*max(abs(uF_true(pick)-U0n(pick)));
end  % iteration loop

uF_true=uF_true
U_fine

L2NormError=L2NormError'
LInfNormError=LInfNormError'
SupperlinearErrorBound=SupperlinearErrorBound'
LinearErrorBound=LinearErrorBound'

% Applying Krylov method
U_old = zeros(size(U_fine));
rhs=UpdateODE(f,U_old);            % your program computes by linearity M^{-1}*f
MinvAfun=@(x) -UpdateODE(U_old,x);  % function that computes M^{-1}*A*x by linearity
% [U_gmres,FLAG,RELRES,ITER,RESVEC] = gmres(MinvAfun,rhs);



[U_gmres,FLAG,RELRES,ITER,RESVEC,L2NormErrorGMRES,LInfNormErrorGMRES] = mygmresODE(MinvAfun,rhs,[],[],[],[],[],U0_fine)
%    [U_gmres, error,RESVEC, iter, flag,L2NormErrorGMRES,LInfNormErrorGMRES] = MygmresODEnew( A, U0_fine, f, M, 10, 10, 1e-12 )

res_Parareal
%Residual of Parareal
figure
% plot(res_Parareal)
% title('Residual of Parareal')

% hold on
%Residual of GMRES
% figure
% plot(RESVEC)
semilogy(0:K,res_Parareal,'*-r',0:length(RESVEC)-1,RESVEC,'*-b','LineWidth',1.5)
legend({'SC two-level additive Schwarz','Parareal with GMRES acceleration'},'Location','northeast','FontSize',12)
xlabel('Iteration')
ylabel('Relative error')
ylim([1e-16 1e2])
% xlim([1e0 1e1+3])
xlim([1e0 20])
title(['T=',num2str(T)])
% title('Residual of Parareal and GMRES')
 set(gca,'xscale','log')
 
% figure
% semilogy(0:K-1,LInfNormError,'b--^',0:K-1,L2NormError,'g--^',0:K-1,SupperlinearErrorBound,'r*-',0:K-1,LinearErrorBound,'mx-')
% legend('L^{\infty}NormError','L^2NormError','Superlinear bound','Linear bound')
% xlabel('k')
% ylim([1e-18 1])
% title('Convergence for Parareal')


% figure
% semilogy(0:length(L2NormErrorGMRES)-1,L2NormErrorGMRES,'g--^',0:length(LInfNormErrorGMRES)-1,LInfNormErrorGMRES,'c--^',0:K-1,LInfNormError,'b--^',0:K-1,L2NormError,'m--^',0:K-1,SupperlinearErrorBound,'r*-',0:K-1,LinearErrorBound,'kx-')
% legend('L^2NormError-GMRES','L^{\infty}NormError-GMRES','L^2NormError-Parareal','L^{\infty}NormError-Parareal','Superlinear bound','Linear bound')
% xlabel('Iteration')
% ylabel('Error')
% ylim([1e-18 1])
% % title(['Convergence for Parareal and GMRES , T=',num2str(T)])
% % set(gca,'xscale','log')



% figure
% semilogy(0:length(L2NormErrorGMRES)-1,L2NormErrorGMRES,'g--^',0:length(LInfNormErrorGMRES)-1,LInfNormErrorGMRES,'c--^',0:K,L2NormError,'m--^',0:K,LInfNormError,'b--^',0:K,SupperlinearErrorBound,'r*-','LineWidth',1.5)
% legend({'L^2- Parareal with GMRES','L^{\infty}- Parareal with GMRES','L^2- FC- Two-level DD','L^{\infty}- FC- Two-level DD','Superlinear bound'},'Location','southwest','FontSize',12) % southwest northeast
% xlabel('Iteration')
% ylabel('Error')
% % ylim([1e-20 1e4])
% ylim([1e-17 1e4])
% % xlim([0 8])
% title(['T=',num2str(T)])
% % title(['Convergence for Parareal and GMRES , T=',num2str(T)])
% % set(gca,'xscale','log')


figure
semilogy(0:length(L2NormErrorGMRES)-1,L2NormErrorGMRES,'g--^',0:length(LInfNormErrorGMRES)-1,LInfNormErrorGMRES,'c--^',0:K,L2NormError,'m--^',0:K,LInfNormError,'b--^','LineWidth',1.5)
legend({'L^2- Parareal with GMRES acceleration','L^{\infty}- Parareal with GMRES acceleration','L^2- SC two-level additive Schwarz','L^{\infty}- SC two-level additive Schwarz'},'Location','northeast','FontSize',12) % southwest northeast
xlabel('Iteration')
ylabel('Error')
% ylim([1e-20 1e4])
ylim([1e-17 1e0])
% xlim([0 8])
title(['T=',num2str(T)])
% title(['Convergence for Parareal and GMRES , T=',num2str(T)])
% set(gca,'xscale','log')

% figure
% semilogy(0:K-1,LInfNormError,'b--^',0:K-1,SupperlinearErrorBound,'r*-',0:K-1,LinearErrorBound,'mx-')
% legend('L^{oo}NormError','Superlinear bound','Linear bound')
% xlabel('k')
% ylim([1e-18 1])
% title(['Convergence for Parareal and GMRES , T=',num2str(T)])

% figure 
% spy(M)
% 
% figure 
% spy(P)
% 
% figure 
% spy(R'*B*R +  Id-R'*R)

RESVEC














