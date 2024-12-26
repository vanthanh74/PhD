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
N = 20;     %   number of coarse time intervals
T = 1;      %   T = N*Dt
Dt = T/N;   %   coarse timesteps

xt = linspace(0,T,N+1);

% fine propagator
m = 2;    % number of fine time steps in each coarse time interval
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
phiDT = (1-a*Dt)^-1;
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

L2NormError = sqrt(sum(Dt*(uF_true(pick)-U_fine(pick)).^2));%norm(uF_true(pick)-U_fine(pick),2);%
LInfNormError = norm(uF_true(pick)-U_fine(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));

SupperlinearErrorBound  = LInfNormError;
LinearErrorBound  = LInfNormError;

% FC_boundL2 = L2NormError;
FC_boundL2 = (abs(phi^m - phiDT)*(1 - abs(phiDT)^((N)/m) )/ (1 - abs(phiDT)))^(0)*sqrt(sum(Dt*(uF_true(pick)-U0n(pick)).^2));%norm(uF_true(pick)-U0n(pick),2);%
   
L2NormErrorParareal = L2NormError;
L2NormErrorNonO2 = L2NormError;

LInfNormErrorParareal = LInfNormError;
LInfNormErrorNonO2= LInfNormError;

FCF_boundL2 = L2NormError;
FCFF_boundL2 = L2NormError;
FCFC_boundL2 = L2NormError;
FCFCF_boundL2 = L2NormError;
%%
% rhs
f = zeros(N_fine_nodes,1);
f(1) = u0;


% phi   = (1-a*dt)^-1;
% phiDT = (1-a*Dt)^-1;

% initial coarse integration Parareal


Uk_0 = zeros(N,1);
M_tilde = eye(N+1,N+1) + (-G_tilde)*diag(ones(N,1),-1);  
rhs = [u0; (F_tilde-G_tilde)*Uk_0];
U_Parareal = M_tilde\rhs;
f1 = zeros(N+1,1);
f1(1) = u0;

% U_Parareal = mdl(u0,T,N,'coarse');

% iteration loop
for k = 1: K
    % k = k + 1

    %----------------------------------------------------------
    % Parareal solution
    
     U_Parareal = U_Parareal+M_tilde\(f1 -  A_tilde*U_Parareal); %use
     
%        for n = 1 : N 
%            % n = n + 1
%            G0 = mdl(U_Parareal(:,n),Dt,1,'coarse');
%            F0 = mdl(U_Parareal(:,n),Dt,m,'fine');
%            G(:,n) = G0(:,end);
%            F(:,n) = F0(:,end);
%        end % end parallel loop
%        
%   
%        % sequential loop - corrector
%        if k < K
%            for n = 1: N
%                Gn =  mdl(U_Parareal(:,n),Dt,1,'coarse');
%                U_Parareal(:, n + 1) = Gn(:,end) - G(:,n) + F(:,n);
%            end % sequential loop
%        end

   %----------------------------------------------------------        
    % 2-level domain decomposition preconditioner Parareal  
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
 
    plot(ttrue(pick),U_fine(pick),'Color',cmap(k,:),'Marker','o');
    
    %
    %        norm(uF_true(pick)-U,'inf')
    proNj(k+1) = proNj(k)*(N-k);
    L2NormErrorNonO2(k+1) = sqrt(sum(Dt*(uF_true(pick)-U_fine(pick)).^2));%norm(uF_true(pick)-U_fine(pick),2);%
    L2NormErrorParareal(k+1) = sqrt(sum(Dt*(uF_true(pick)-U_Parareal).^2));%norm(uF_true(pick)-U_fine(pick),2);%
    
    LInfNormErrorNonO2(k+1) = norm(uF_true(pick)-U_fine(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));%
    LInfNormErrorParareal(k+1) = norm(uF_true(pick)-U_Parareal,inf);%max(abs(uF_true(pick)-U_fine(pick)));%
   
    SupperlinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt)))^(k)/factorial(k)*proNj(k+1)*max(abs(uF_true(pick)-U0n(pick)));
    LinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt))/(1- abs(Rz(a*Dt))))^(k)*max(abs(uF_true(pick)-U0n(pick)));
    
    %FC bound
    FC_boundL2(k+1) = (abs(phi^m - phiDT)*(1 - abs(phiDT)^((N)/m) )/ (1 - abs(phiDT)))^(k)*sqrt(sum(Dt*(uF_true(pick)-U0n(pick)).^2));%norm(uF_true(pick)-U0n(pick),2);%
   
    FCF_boundL2(k+1) = (abs(phi^m - phiDT)*(1 - abs(phiDT)^((N)/m-1) )/ (1 - abs(phiDT))*abs(phi)^(m))^(k)*sqrt(sum(Dt*(uF_true(pick)-U0n(pick)).^2));
    FCFF_boundL2(k+1) = (abs(phi^m - phiDT)*(1 - abs(phiDT)^((N)/m-2) )/ (1 - abs(phiDT))*abs(phi)^(2*m))^(k)*sqrt(sum(Dt*(uF_true(pick)-U0n(pick)).^2));
%     FCFC_boundL2(k+1) = ((phi^m - phiDT)^2*(1 - ((N)/m )*abs(phiDT)^((N)/m-1) + ((N)/m-1)*abs(phiDT)^((N-1)/m) )/ (1 - abs(phiDT))^2)^(k)*sqrt(sum(Dt*(uF_true(pick)-U0n(pick)).^2));
    FCFCF_boundL2(k+1) = ((phi^m - phiDT)^2*(1 - ((N)/m - 1 )*abs(phiDT)^((N)/m-2) + ((N)/m-2)*abs(phiDT)^((N-1)/m-1) )/ (1 - abs(phiDT))^2*abs(phi)^m)^(k)*sqrt(sum(Dt*(uF_true(pick)-U0n(pick)).^2));

end  % iteration loop

uF_true=uF_true
U_fine

% L2NormError=L2NormError'
% LInfNormError=LInfNormError'
% SupperlinearErrorBound=SupperlinearErrorBound'
% LinearErrorBound=LinearErrorBound'
% FC_boundL2 = FC_boundL2'
% FCF_boundL2 = FCF_boundL2'
% FCFF_boundL2 = FCFF_boundL2'
% FCFC_boundL2 = FCFC_boundL2'
% FCFCF_boundL2 = FCFCF_boundL2'


figure
semilogy(0:K,LInfNormErrorParareal,'b--^',0:K,L2NormErrorParareal,'g--^',0:K,SupperlinearErrorBound,'r*-',0:K,LinearErrorBound,'mx-',0:K,FC_boundL2,'k.-',0:K,FCF_boundL2,'y.-')
legend('L^{\infty}NormError','L^2NormError','Superlinear bound','Linear bound','FC boundL2','FCF boundL2')
xlabel('k')
ylim([1e-18 1])

%----------------------------
% %FC bound L2
% % % superlinear with Gander's bound
% figure
% semilogy(0:K,LInfNormErrorParareal,'g--^',0:K,L2NormErrorParareal,'m--^',0:K,L2NormErrorNonO2,'r--s',0:K,LInfNormErrorNonO2,'b--s',0:K,SupperlinearErrorBound,'ko-',0:K,FC_boundL2,'co-','LineWidth',1.5)
% legend({'L^{\infty}- Parareal','L^2- Parareal','L^2- SC two-level additive Schwarz','L^{\infty}- SC two-level additive Schwarz','Superlinear bound','SC bound L^2- norm'},'Location','northeast','FontSize',12)  %southwest % northeast
% xlabel('Iteration ')
% ylabel('Error')
% % ylim([1e-18 1e-1])
% % ylim([1e-14 1e4])
% ylim([1e-22 1e0])
% % xlim([0 7])
% title(['T=',num2str(T)])

%%%% superlinear without Gander's bound
% figure
% semilogy(0:K,LInfNormErrorParareal,'g--^',0:K,L2NormErrorParareal,'m--^',0:K,L2NormErrorNonO2,'r--s',0:K,LInfNormErrorNonO2,'b--s',0:K,FC_boundL2,'co-','LineWidth',1.5)
% legend({'L^{\infty}- Parareal','L^2- Parareal','L^2- SC two-level additive Schwarz','L^{\infty}- SC two-level additive Schwarz','SC bound L^2- norm'},'Location','northeast','FontSize',12)  %southwest % northeast
% xlabel('Iteration ')
% ylabel('Error')
% % ylim([1e-18 1e-1])
% % ylim([1e-14 1e4])
% ylim([1e-18 1e-1])
% xlim([0 7])
% title(['T=',num2str(T)])


%%%% superlinear without Gander's bound and 2-norm
figure
semilogy(0:K,L2NormErrorParareal,'b--^',0:K,L2NormErrorNonO2,'r--s',0:K,FC_boundL2,'ko-','LineWidth',1.5)
legend({'Parareal','SC two-level additive Schwarz','SC bound'},'Location','northeast','FontSize',12)  %southwest % northeast
xlabel('Iteration ')
ylabel('Error')
% ylim([1e-18 1e-1])
% ylim([1e-14 1e4])
ylim([1e-18 1e-1])
xlim([0 7])
title(['T=',num2str(T)])


% % % linear with Gander's bound
% figure
% semilogy(0:K,LInfNormErrorParareal,'g--^',0:K,L2NormErrorParareal,'m--^',0:K,L2NormErrorNonO2,'r--s',0:K,LInfNormErrorNonO2,'b--s',0:K,LinearErrorBound,'ko-',0:K,FC_boundL2,'co-','LineWidth',1.5)
% legend({'L^{\infty}- Parareal','L^2- Parareal','L^2- SC two-level additive Schwarz','L^{\infty}- SC two-level additive Schwarz','Linear bound','SC bound L^2- norm'},'Location','southwest','FontSize',12)  %southwest % northeast
% xlabel('Iteration ')
% ylabel('Error')
% % ylim([1e-18 1e-1])
% % ylim([1e-14 1e4])
% ylim([1e-16 1e0])
% % xlim([0 7])
% title(['T=',num2str(T)])
% 
% % % % linear without Gander's bound
% figure
% semilogy(0:K,LInfNormErrorParareal,'g--^',0:K,L2NormErrorParareal,'m--^',0:K,L2NormErrorNonO2,'r--s',0:K,LInfNormErrorNonO2,'b--s',0:K,LinearErrorBound,'ko-',0:K,FC_boundL2,'co-','LineWidth',1.5)
% legend({'L^{\infty}- Parareal','L^2- Parareal','L^2- SC two-level additive Schwarz','L^{\infty}- SC two-level additive Schwarz','Linear bound','SC bound L^2- norm'},'Location','southwest','FontSize',12)  %southwest % northeast
% xlabel('Iteration ')
% ylabel('Error')
% % ylim([1e-18 1e-1])
% % ylim([1e-14 1e4])
% ylim([1e-16 1e0])
% % xlim([0 7])
% title(['T=',num2str(T)])

