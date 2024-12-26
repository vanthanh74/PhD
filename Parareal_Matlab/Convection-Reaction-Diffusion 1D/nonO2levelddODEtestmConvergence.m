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
% L2NormError = norm(uF_true(pick)-U_fine(pick),2);
LInfNormError = norm(uF_true(pick)-U_fine(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));

SupperlinearErrorBound  = LInfNormError;
LinearErrorBound  = LInfNormError;
FC_boundL2 = L2NormError;
FCF_boundL2 = L2NormError;
FCFF_boundL2 = L2NormError;
FCFC_boundL2 = L2NormError;
FCFCF_boundL2 = L2NormError;
%%
% rhs
f = zeros(N_fine_nodes,1);
f(1) = u0;


phi   = (1-a*dt)^-1;
phiDT = (1-a*Dt)^-1;
E = [];
% iteration loop
for k = 1: K-1
    % k = k + 1
    %         rhs=[u0];
    %         for i=1:2:2*N
    %             rhs=[rhs; phi*U_fine(i);-G_tilde*U_fine(i)];
    %         end
    %         U_fine = M\rhs;
    %         U_fine = U_fine + M^(-1)*(f-A*U_fine);
    
%     r_k = f - A*U_fine;
%     U_fine_temp = P\r_k;
%     U_fine = U_fine + (R'*B*R +  Id-R'*R)\U_fine_temp;

    
%     r_k = f - A*U_fine;
%     U_fine_temp = P\r_k;
%     U_fine = U_fine + (R'*B*R +  Id-R'*R)\U_fine_temp;
    
%     U_fine = U_fine + (R'*B*R +  Id-R'*R)^(-1)*P^(-1)*r_k;
    
%     U_fine = M\(M*U_fine + (f-A*U_fine)); % use
    
    
%      U_fine = M\(M*U_fine + (f-A*U_fine));
    %         U_fine = gmres(M,M*U_fine + (f-A*U_fine),10,1e-6);
    
    
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
    
   
    %-------------------------------------------------
% %     % FCF 
% % %     % residual
%     r_k = f - A*U_fine;
% % %     first level ~ fine propagator ~ additive Schwarz in parallel
% % %         U_fine_temp = P\r_k;   
%     for i = 1:N+1
%         if i == 1
%             x_temp_1 = R1*r_k;
%             x_temp_2 = A1\x_temp_1;
%             x{i} = R1'*x_temp_2;
%             U_fine_temp1 = x{i};         
%         else
%             x_temp_1 = R_sub{i}*r_k;
%             x_temp_2 = A_sub{i}\x_temp_1;
%             x{i} = R_sub{i}'*x_temp_2;  
%             U_fine_temp1 = U_fine_temp1 + x{i};
%         end            
%     end
%     
%         U_fine = U_fine + U_fine_temp1;
%         
% % %     -------------------------------------------------
% %      %FCFF   
% % % %         residual
%     r_k = f - A*U_fine;
% % %     first level ~ fine propagator ~ additive Schwarz in parallel
% % %         U_fine_temp = P\r_k;   
%     for i = 1:N+1
%         if i == 1
%             x_temp_1 = R1*r_k;
%             x_temp_2 = A1\x_temp_1;
%             x{i} = R1'*x_temp_2;
%             U_fine_temp2 = x{i};         
%         else
%             x_temp_1 = R_sub{i}*r_k;
%             x_temp_2 = A_sub{i}\x_temp_1;
%             x{i} = R_sub{i}'*x_temp_2;  
%             U_fine_temp2 = U_fine_temp2 + x{i};
%         end            
%     end
%     
%         U_fine = U_fine + U_fine_temp2;

    %-------------------------------------------------  
% % %      %FCFC
%             % residual
%     r_k = f - A*U_fine;
%     
% %   % first level ~ fine propagator ~ additive Schwarz in parallel
% % %     U_fine_temp = P\r_k;   
%     for i = 1:N+1
% % %         i = i + 1
%         if i == 1
%             x_temp_1 = R1*r_k;
%             x_temp_2 = A1\x_temp_1;
%             x{i} = R1'*x_temp_2;
%             U_fine_temp = x{i};         
%         else
%             x_temp_1 = R_sub{i}*r_k;
%             x_temp_2 = A_sub{i}\x_temp_1;
%             x{i} = R_sub{i}'*x_temp_2;  
%             U_fine_temp = U_fine_temp + x{i};
%         end            
%     end
% 
% % %   % second level ~ coarse grid correction sequentially
%     U_fine_temp_1 = R*U_fine_temp;
%     U_fine_temp_2 = B\U_fine_temp_1;
%     U_fine_temp_3 = R'*U_fine_temp_2;
%     U_fine_temp_4 = U_fine_temp_3 + U_fine_temp - R'*R*U_fine_temp;
%     U_fine = U_fine + U_fine_temp_4;
%     
     %-------------------------------------------------  
% % %      %FCFCF
% % %         % residual
%     r_k = f - A*U_fine;
% % %     % first level ~ fine propagator ~ additive Schwarz in parallel
% % %     %     U_fine_temp = P\r_k;   
%     for i = 1:N+1
%         if i == 1
%             x_temp_1 = R1*r_k;
%             x_temp_2 = A1\x_temp_1;
%             x{i} = R1'*x_temp_2;
%             U_fine_temp3 = x{i};         
%         else
%             x_temp_1 = R_sub{i}*r_k;
%             x_temp_2 = A_sub{i}\x_temp_1;
%             x{i} = R_sub{i}'*x_temp_2;  
%             U_fine_temp3 = U_fine_temp3 + x{i};
%         end            
%     end
%     
%         U_fine = U_fine + U_fine_temp3;
    
    
    
    
    
    
    
    
    
    
    
    
    plot(ttrue(pick),U_fine(pick),'Color',cmap(k,:),'Marker','o');
    
    %
    %        norm(uF_true(pick)-U,'inf')
    proNj(k+1) = proNj(k)*(N-k);
    L2NormError(k+1) = sqrt(sum(Dt*(uF_true(pick)-U_fine(pick)).^2));%norm(uF_true(pick)-U_fine(pick),2);%
%     L2NormError(k+1) = norm(uF_true(pick)-U_fine(pick),2);%max(abs(uF_true(pick)-U_fine(pick)));%   
    LInfNormError(k+1) = norm(uF_true(pick)-U_fine(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));%
    SupperlinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt)))^(k)/factorial(k)*proNj(k+1)*max(abs(uF_true(pick)-U0n(pick)));
    LinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt))/(1- abs(Rz(a*Dt))))^(k)*max(abs(uF_true(pick)-U0n(pick)));
    
    %FC bound
%     FC_boundL2(k+1) = (abs(phi^m - phiDT)*(1 - abs(phiDT)^((N-1)/m) )/ (1 - abs(phiDT)))^(k)*norm(abs(uF_true(pick)-U0n(pick)),2);
      FC_boundL2(k+1) = (abs(phi^m - phiDT)*(1 - abs(phiDT)^((N-1)/m) )/ (1 - abs(phiDT)))^(k)*sqrt(sum(Dt*(uF_true(pick)-U0n(pick)).^2));

    FCF_boundL2(k+1) = (abs(phi^m - phiDT)*(1 - abs(phiDT)^((N-1)/m-1) )/ (1 - abs(phiDT))*abs(phi)^(m))^(k)*sqrt(sum(Dt*(uF_true(pick)-U0n(pick)).^2));
    FCFF_boundL2(k+1) = (abs(phi^m - phiDT)*(1 - abs(phiDT)^((N-1)/m-2) )/ (1 - abs(phiDT))*abs(phi)^(2*m))^(k)*sqrt(sum(Dt*(uF_true(pick)-U0n(pick)).^2));
    FCFC_boundL2(k+1) = ((phi^m - phiDT)^2*(1 - ((N-1)/m )*abs(phiDT)^((N-1)/m-1) + ((N-1)/m-1)*abs(phiDT)^((N-1)/m) )/ (1 - abs(phiDT))^2)^(k)*sqrt(sum(Dt*(uF_true(pick)-U0n(pick)).^2));
    FCFCF_boundL2(k+1) = ((phi^m - phiDT)^2*(1 - ((N-1)/m - 1 )*abs(phiDT)^((N-1)/m-2) + ((N-1)/m-2)*abs(phiDT)^((N-1)/m-1) )/ (1 - abs(phiDT))^2*sqrt(sum(Dt*(uF_true(pick)-U0n(pick)).^2)));

%       FCF_boundL2(k+1) = (abs(phi^m - phiDT)*(1 - abs(phiDT)^((N-1)/m-1) )/ (1 - abs(phiDT))*abs(phi)^(m))^(k)*norm(abs(uF_true(pick)-U0n(pick)),2);
%     FCFF_boundL2(k+1) = (abs(phi^m - phiDT)*(1 - abs(phiDT)^((N-1)/m-2) )/ (1 - abs(phiDT))*abs(phi)^(2*m))^(k)*norm(abs(uF_true(pick)-U0n(pick)),2);
%     FCFC_boundL2(k+1) = ((phi^m - phiDT)^2*(1 - ((N-1)/m )*abs(phiDT)^((N-1)/m-1) + ((N-1)/m-1)*abs(phiDT)^((N-1)/m) )/ (1 - abs(phiDT))^2)^(k)*norm(abs(uF_true(pick)-U0n(pick)),2);
%     FCFCF_boundL2(k+1) = ((phi^m - phiDT)^2*(1 - ((N-1)/m - 1 )*abs(phiDT)^((N-1)/m-2) + ((N-1)/m-2)*abs(phiDT)^((N-1)/m-1) )/ (1 - abs(phiDT))^2*norm(abs(uF_true(pick)-U0n(pick)),2));
E = [E; (abs(phi^m - phiDT)*(1 - abs(phiDT)^((N-1)/m) )/ (1 - abs(phiDT)))^(k)]

end  % iteration loop
phi = phi
phiDT = phiDT
phi_m = phi^m
% uF_true=uF_true
% U_fine
% 
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
semilogy(0:K-1,LInfNormError,'b--^',0:K-1,L2NormError,'g--^',0:K-1,SupperlinearErrorBound,'r*-',0:K-1,LinearErrorBound,'mx-',0:K-1,FC_boundL2,'k.-',0:K-1,FCF_boundL2,'y.-')
legend('L^{\infty}NormError','L^2NormError','Superlinear bound','Linear bound','FC boundL2','FCF boundL2')
xlabel('k')
ylim([1e-18 1])

%----------------------------
% %FC bound L2
figure
semilogy(0:K-1,L2NormError,'b--^',0:K-1,SupperlinearErrorBound,'k.-',0:K-1,FC_boundL2,'ro-')
legend('L^2- 2-level DD error norm','Superlinear bound','FC bound L^2- norm')
xlabel('Iteration ')
ylabel('Error')
ylim([1e-18 1])
xlim([0 7])
title(['T=',num2str(T)])

figure
semilogy(0:K-1,LInfNormError,'m--^',0:K-1,L2NormError,'b--^',0:K-1,SupperlinearErrorBound,'k.-',0:K-1,FC_boundL2,'ro-')
legend('L^{\infty}NormError','L^2- 2-level DD error norm','Superlinear bound L^{\infty}- norm','FC bound L^2- norm')
xlabel('Iteration ')
ylabel('Error')
% ylim([1e-18 1])
% xlim([0 7])
title(['T=',num2str(T)])
%----------------------------
%FCF bound L2
% figure
% semilogy(0:K-1,L2NormError,'b--^',0:K-1,SupperlinearErrorBound,'k.-',0:K-1,FCF_boundL2,'ro-')
% legend('L^2NormError','Superlinear bound','FCF boundL2')
% xlabel('k')
% ylim([1e-18 1])
% title(['T=',num2str(T)])
%----------------------------
%FCFF bound L2
% figure
% semilogy(0:K-1,L2NormError,'b--^',0:K-1,SupperlinearErrorBound,'k.-',0:K-1,FCFF_boundL2,'ro-')
% legend('L^2NormError','Superlinear bound','FCFF boundL2')
% xlabel('k')
% ylim([1e-18 1])
% title(['T=',num2str(T)])



%----------------------------
% %FCFC bound L2
% figure
% semilogy(0:K-1,L2NormError,'b--^',0:K-1,SupperlinearErrorBound,'k.-',0:K-1,FCFC_boundL2,'ro-')
% legend('L^2NormError','Superlinear bound','FCFC boundL2')
% xlabel('k')
% ylim([1e-18 1])
% title(['T=',num2str(T)])

% %----------------------------
% % %FCFC bound L2
% figure
% semilogy(0:K-1,L2NormError,'b--^',0:K-1,SupperlinearErrorBound,'k.-',0:K-1,FCFCF_boundL2,'ro-')
% legend('L^2NormError','Superlinear bound','FCFCF boundL2')
% xlabel('k')
% ylim([1e-18 1])
% title(['T=',num2str(T)])






% figure
% semilogy(0:K-1,LInfNormError,'b--^',0:K-1,SupperlinearErrorBound,'r*-',0:K-1,LinearErrorBound,'mx-')
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
% 2*(phi^m-phiDT)*phiDT*(phi^m-phiDT)-2*phiDT*(phi^m-phiDT)^2



