close all
clear
clc


%% ---------------------------------------------------------------
%  Parareal algorithm for linear problem:
%         y'(t) = a*y(t), t in (0,T)
%          y(t=0) = y0


% non-overlapping subdomains 2-level domain decomposition  preconditioner

global a ;
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
dt = Dt/m

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

figure
hold on
plot(ttrue,uF_true,'-^k',ttrue,U_fine,'-*m');
legend('True solution','Parareal solution at k = 0')

U0n = U_fine;
proNj = 1;

L2NormError_FC = sqrt(sum(Dt*(uF_true(pick)-U_fine(pick)).^2));%norm(uF_true(pick)-U_fine(pick),2);%
LInfNormError_FC = norm(uF_true(pick)-U_fine(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));


L2NormError_FCF = sqrt(sum(Dt*(uF_true(pick)-U_fine(pick)).^2));%norm(uF_true(pick)-U_fine(pick),2);%
LInfNormError_FCF = norm(uF_true(pick)-U_fine(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));

L2NormError_FCFF = sqrt(sum(Dt*(uF_true(pick)-U_fine(pick)).^2));%norm(uF_true(pick)-U_fine(pick),2);%
LInfNormError_FCFF = norm(uF_true(pick)-U_fine(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));

L2NormError_FCFC = sqrt(sum(Dt*(uF_true(pick)-U_fine(pick)).^2));%norm(uF_true(pick)-U_fine(pick),2);%
LInfNormError_FCFC = norm(uF_true(pick)-U_fine(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));

L2NormError_FCFCF = sqrt(sum(Dt*(uF_true(pick)-U_fine(pick)).^2));%norm(uF_true(pick)-U_fine(pick),2);%
LInfNormError_FCFCF = norm(uF_true(pick)-U_fine(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));

SupperlinearErrorBound  = L2NormError_FC;
LinearErrorBound  = L2NormError_FC;

FC_boundL2 = L2NormError_FC;
FCF_boundL2 = L2NormError_FC;
FCFF_boundL2 = L2NormError_FC;
FCFC_boundL2 = L2NormError_FC;
FCFCF_boundL2 = L2NormError_FC;


%%
% rhs
f = zeros(N_fine_nodes,1);
f(1) = u0;


phi   = (1-a*dt)^-1;
phiDT = (1-a*Dt)^-1;

U_FC = U_fine;
U_FCF = U_fine;
U_FCFF = U_fine;
U_FCFC = U_fine;
U_FCFCF = U_fine;

% iteration loop
for k = 1: K
    % k = k + 1
 
    % FC - relaxation
    [U_FC,L2NormError_FC_k,LInfNormError_FC_k] = TwolevelDDsolveODE(f,A,U_FC,N,R1,A1,R_sub,A_sub,R,B,pick,Dt,uF_true,'FC')
    L2NormError_FC(k+1) = L2NormError_FC_k;
    LInfNormError_FC(k+1) = LInfNormError_FC_k;
    
      % FCF - relaxation
    [U_FCF,L2NormError_FCF_k,LInfNormError_FCF_k] = TwolevelDDsolveODE(f,A,U_FCF,N,R1,A1,R_sub,A_sub,R,B,pick,Dt,uF_true,'FCF')
    L2NormError_FCF(k+1) = L2NormError_FCF_k;
    LInfNormError_FCF(k+1) = LInfNormError_FCF_k;
    
          % FCFF - relaxation
    [U_FCFF,L2NormError_FCFF_k,LInfNormError_FCFF_k] = TwolevelDDsolveODE(f,A,U_FCFF,N,R1,A1,R_sub,A_sub,R,B,pick,Dt,uF_true,'FCFF')
    L2NormError_FCFF(k+1) = L2NormError_FCFF_k;
    LInfNormError_FCFF(k+1) = LInfNormError_FCFF_k;
    
    
             % FCFC - relaxation
    [U_FCFC,L2NormError_FCFC_k,LInfNormError_FCFC_k] = TwolevelDDsolveODE(f,A,U_FCFC,N,R1,A1,R_sub,A_sub,R,B,pick,Dt,uF_true,'FCFC')
    L2NormError_FCFC(k+1) = L2NormError_FCFC_k;
    LInfNormError_FCFC(k+1) = LInfNormError_FCFC_k;
    
                 % FCFCF - relaxation
    [U_FCFCF,L2NormError_FCFCF_k,LInfNormError_FCFCF_k] = TwolevelDDsolveODE(f,A,U_FCFCF,N,R1,A1,R_sub,A_sub,R,B,pick,Dt,uF_true,'FCFCF')
    L2NormError_FCFCF(k+1) = L2NormError_FCFCF_k;
    LInfNormError_FCFCF(k+1) = LInfNormError_FCFCF_k;
    
    
%     plot(ttrue(pick),U_fine(pick),'Color',cmap(k,:),'Marker','o');
    
    %
    %        norm(uF_true(pick)-U,'inf')
    proNj(k+1) = proNj(k)*(N-k);
%     L2NormError(k+1) = sqrt(sum(Dt*(uF_true(pick)-U_fine(pick)).^2));%norm(uF_true(pick)-U_fine(pick),2);%
%     LInfNormError(k+1) = norm(uF_true(pick)-U_fine(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));%
    SupperlinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt)))^(k)/factorial(k)*proNj(k+1)*max(abs(uF_true(pick)-U0n(pick)));
    LinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt))/(1- abs(Rz(a*Dt))))^(k)*max(abs(uF_true(pick)-U0n(pick)));
    
    %FC bound
    FC_boundL2(k+1) = (abs(phi^m - phiDT)*(1 - abs(phiDT)^((N)/m) )/ (1 - abs(phiDT)))^(k)*sqrt(sum(Dt*(uF_true(pick)-U0n(pick)).^2));
    FCF_boundL2(k+1) = (abs(phi^m - phiDT)*(1 - abs(phiDT)^((N)/m-1) )/ (1 - abs(phiDT))*abs(phi)^(m))^(k)*sqrt(sum(Dt*(uF_true(pick)-U0n(pick)).^2));
    FCFF_boundL2(k+1) = (abs(phi^m - phiDT)*(1 - abs(phiDT)^((N)/m-2) )/ (1 - abs(phiDT))*abs(phi)^(2*m))^(k)*sqrt(sum(Dt*(uF_true(pick)-U0n(pick)).^2));
    FCFC_boundL2(k+1) = ((phi^m - phiDT)^2*(1 - ((N)/m )*abs(phiDT)^((N)/m-1) + ((N)/m-1)*abs(phiDT)^((N)/m) )/ (1 - abs(phiDT))^2)^(k)*sqrt(sum(Dt*(uF_true(pick)-U0n(pick)).^2));
    FCFCF_boundL2(k+1) = ((phi^m - phiDT)^2*(1 - ((N)/m - 1 )*abs(phiDT)^((N)/m-2) + ((N)/m-2)*abs(phiDT)^((N)/m-1) )/ (1 - abs(phiDT))^2*abs(phi)^m)^(k)*sqrt(sum(Dt*(uF_true(pick)-U0n(pick)).^2));

end  % iteration loop

uF_true=uF_true
U_fine

L2NormError_FC=L2NormError_FC'
LInfNormError_FC=LInfNormError_FC'

L2NormError_FCF=L2NormError_FCF'
LInfNormError_FCF=LInfNormError_FCF'

L2NormError_FCFF=L2NormError_FCFF'
LInfNormError_FCFF=LInfNormError_FCFF'

L2NormError_FCFC=L2NormError_FCFC'
LInfNormError_FCFC=LInfNormError_FCFC'

L2NormError_FCFCF=L2NormError_FCFCF'
LInfNormError_FCFCF=LInfNormError_FCFCF'

SupperlinearErrorBound=SupperlinearErrorBound'
LinearErrorBound=LinearErrorBound'
FC_boundL2 = FC_boundL2'
FCF_boundL2 = FCF_boundL2'
FCFF_boundL2 = FCFF_boundL2'
FCFC_boundL2 = FCFC_boundL2'
FCFCF_boundL2 = FCFCF_boundL2'


% figure
% semilogy(0:K-1,LInfNormError,'b--^',0:K-1,L2NormError,'g--^',0:K-1,SupperlinearErrorBound,'r*-',0:K-1,LinearErrorBound,'mx-',0:K-1,FC_boundL2,'k.-',0:K-1,FCF_boundL2,'y.-')
% legend('L^{\infty}NormError','L^2NormError','Superlinear bound','Linear bound','FC boundL2','FCF boundL2')
% xlabel('k')
% ylim([1e-18 1])

%----------------------------
% %FC bound L2
figure
semilogy(0:K,L2NormError_FC,'b--^',0:K,SupperlinearErrorBound,'k.-',0:K,FC_boundL2,'ro-')
legend('L^2NormError','Superlinear bound','FC boundL2')
xlabel('k')
ylim([1e-18 1])
title(['T=',num2str(T)])
%----------------------------
% %FCF bound L2
figure
semilogy(0:K,L2NormError_FCF,'b--^',0:K,SupperlinearErrorBound,'k.-',0:K,FCF_boundL2,'ro-')
legend('L^2NormError','Superlinear bound','FCF boundL2')
xlabel('k')
ylim([1e-18 1])
title(['T=',num2str(T)])
%----------------------------
%FCFF bound L2
figure
semilogy(0:K,L2NormError_FCFF,'b--^',0:K,SupperlinearErrorBound,'k.-',0:K,FCFF_boundL2,'ro-')
legend('L^2NormError','Superlinear bound','FCFF boundL2')
xlabel('k')
ylim([1e-18 1])
title(['T=',num2str(T)])



%----------------------------
% %FCFC bound L2
figure
semilogy(0:K,L2NormError_FCFC,'b--^',0:K,SupperlinearErrorBound,'k.-',0:K,FCFC_boundL2,'ro-')
legend('L^2NormError','Superlinear bound','FCFC boundL2')
xlabel('Iteration')
ylabel('Error')
ylim([1e-18 1])
title(['T=',num2str(T)])

% %----------------------------
% % %FCFC bound L2
figure
semilogy(0:K,L2NormError_FCFCF,'b--^',0:K,SupperlinearErrorBound,'k.-',0:K,FCFCF_boundL2,'ro-')
legend('L^2NormError','Superlinear bound','FCFCF boundL2')
xlabel('Iteration')
ylabel('Error')
ylim([1e-18 1])
title(['T=',num2str(T)])

% %----------------------------
% figure
% semilogy(0:K,L2NormError_FC,'b--^',0:K,L2NormError_FCF,'m--^',0:K,L2NormError_FCFF,'c--^',0:K,L2NormError_FCFC,'g--^',0:K,L2NormError_FCFCF,'r--^','LineWidth',1.5)
% legend({'L^2- FC- Two-level DD','L^2- FCF- Two-level DD','L^2- FCFF- Two-level DD','L^2- FCFC- Two-level DD','L^2- FCFCF- Two-level DD'},'Location','northeast','FontSize',12)
% xlabel('Iteration')
% ylabel('Error')
% % ylim([1e-18 1e-1])
% title(['T=',num2str(T)])

figure
semilogy(0:K,L2NormError_FC,'b--^',0:K,L2NormError_FCF,'m--^',0:K,L2NormError_FCFF,'c--^',0:K,L2NormError_FCFCF,'r--^','LineWidth',1.5) % southwest
legend({'SC two-level additive Schwarz','SCS two-level additive Schwarz','SCS^{2} two-level additive Schwarz','S(CS)^{2} two-level additive Schwarz'},'Location','southwest','FontSize',12) % northeast southwest
xlabel('Iteration')
ylabel('Error')
% ylim([1e-18 1e-1])
title(['T=',num2str(T)])





% figure
% semilogy(0:K,L2NormError_FC,'b-^',0:K,L2NormError_FCF,'m-^',0:K,L2NormError_FCFF,'c-^',0:K,L2NormError_FCFC,'g-^',0:K,L2NormError_FCFCF,'r-^')
% legend('L^2NormError FC-relaxation','L^2NormError FCF-relaxation','L^2NormError FCFF-relaxation','L^2NormError FCFC-relaxation','L^2NormError FCFCF-relaxation')
% xlabel('k')
% ylim([1e-18 1])
% title(['T=',num2str(T)])

% figure
% semilogy(0:K,LInfNormError,'b--^',0:K,SupperlinearErrorBound,'r*-',0:K,LinearErrorBound,'mx-')
% legend('L^{oo}NormError','Superlinear bound','Linear bound')
% xlabel('k')
% ylim([1e-18 1])
% title(['T=',num2str(T)])

% figure 
% spy(M)
% 
% figure 
% spy(P)
% 
% figure 
% spy(R'*B*R +  Id-R'*R)

%F-relaxation
% F_rex_error_propagation = eye(2*N+1)-M^-1*A



% phi   = (1-a*dt)^-1;
% phiDT = (1-a*Dt)^-1;
% eig_F = eig(phi)
% eig_C = eig(phiDT)
% 
% L2NormError_FC = L2NormError
% FC_bound = max(abs(phi^m - phiDT)*(1 - abs(phiDT)^((N-1)/m) )/ (1 - abs(phiDT)))*max(abs(uF_true(pick)-U0n(pick)))
% 
% FCF_bound = max(abs(phi^m - phiDT)*(1 - abs(phiDT)^((N-1)/m) )/ (1 - abs(phiDT))*abs(phi)^m)
% 
% SupperlinearErrorBound
% 
% LinearErrorBound



% test commute
2*(phi^m-phiDT)*phiDT*(phi^m-phiDT)-2*phiDT*(phi^m-phiDT)^2



