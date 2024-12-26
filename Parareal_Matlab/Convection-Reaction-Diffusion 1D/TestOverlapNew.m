close all
clear
clc


%% ---------------------------------------------------------------
%  Parareal algorithm for linear problem:  
%         y'(t) = a*y(t), t in (0,T)
%          y(t=0) = y0


% overlapping subdomains 2-level domain decomposition additive preconditioner

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
N = 2;     %   number of coarse time intervals
T = 1;      %   T = N*Dt  
Dt = T/N;   %   coarse timesteps

xt = linspace(0,T,N+1);

% fine propagator
m = 2   % number of fine time steps in each coarse time interval
dt = Dt/m;

K = 7 % number of Parareal iterations

cmap = hsv(K); % color map

% 2-level domain decomposition initial
N_fine_nodes = m*N+1; % 11
N_subdomains = N; % 6
Id = eye(N_fine_nodes);

% plot the time domain
figure 
plot(0:N_fine_nodes-1,1,'.b',0:m:N_fine_nodes-1,1,'ro')
xlim([0 N_fine_nodes-1])
ylim([0.5 1.5])

%%
% original problem matrix
phi = (1-a*dt)^-1;
A = eye(m*N+1,m*N+1) - (phi)*diag(ones(m*N,1),-1);  

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

 xtrue = linspace(0,T,N*m+1);
 pick = 1:m:N*m+1;
 
 uF_true = A\rhs0;
 
% reduced matrix
A_tilde = eye(N+1,N+1) + (-F_tilde)*diag(ones(N,1),-1);      

% approximate reduced matrix
B = eye(N+1,N+1) + (-G_tilde)*diag(ones(N,1),-1);     


% % coarse space initial
indice_R_coarse  = 0:m:m*N;
R = zeros(length(indice_R_coarse),N_fine_nodes);
for i = 1:length(indice_R_coarse)
    R(i,indice_R_coarse(i)+1) = 1;
end



% define nodes in each subdomains
N_nodes_subdomain = 2;

% define overlap region
N_nodes_overlap = 1; % 1 % m/2+1 % m

%define the number of overlapping subdomains
if N_nodes_overlap == 1
    N_subdomains_overlap = N_fine_nodes-1;
elseif N_nodes_overlap == m/2+1
    N_subdomains_overlap = 2*N-1;
elseif N_nodes_overlap == m
    N_subdomains_overlap = N_fine_nodes - N_nodes_overlap;
end

% first subdomain  Omega_1 = {0}
R1 = zeros(1,N_fine_nodes);
R1(1) = 1;
A1 = R1*A*R1';
D1 = 1;


indice_R_subdomains = zeros(N_subdomains_overlap-1,N_nodes_subdomain );
for i =1:N_subdomains_overlap-1 
    indice_R_subdomains(i,:) =  (i-1)*(N_nodes_subdomain-N_nodes_overlap)+N_nodes_overlap:(i-1)*(N_nodes_subdomain-N_nodes_overlap)+N_nodes_overlap+1;
end

% % Omega_1 = {0,1,2} Omega_3 = {2,3,4} 
R_sub_mat = zeros((N_subdomains_overlap-1)*N_nodes_subdomain,N_fine_nodes) ;
for i = 1:N_subdomains_overlap-1
    for j = 1:N_nodes_subdomain
        R_sub_mat(N_nodes_subdomain*(i-1)+j,indice_R_subdomains(i)+j) = 1;
    end
end

R_sub = [];       % restriction matrices
A_sub = [];       % sub-matrices
for i = 2: N_subdomains_overlap
    R_sub{i} = R_sub_mat(N_nodes_subdomain*(i-2)+1:N_nodes_subdomain*(i-2)+N_nodes_subdomain,:);
    A_sub{i} = R_sub{i}*A*R_sub{i}';    
end

% partition of unity
D_sub = [];
for i = 2:N_subdomains_overlap
    % i = i + 1
    for j = 1:N_nodes_subdomain
         D_sub{i}(j,j) = 1/length(find(indice_R_subdomains==indice_R_subdomains(i-1,j)));
    end
end

% % define 2-level domain decomposition additive preconditioner
P = R1'*D1*A1*R1;
D = R1'*D1*R1;
for i = 2:N_subdomains_overlap
    P = P + R_sub{i}'*D_sub{i}*A_sub{i}*R_sub{i};
    D = D + R_sub{i}'*D_sub{i}*R_sub{i};
end
% P = R1'*D1*A1*R1 + R2'*D2*A2*R2 ;
M = P + (R'*B*R - R'*R )
D;







%% ------------------------------------------------------------------

% initial coarse integration

B1 =  eye(N+1,N+1) + (-G_tilde)*diag(ones(N,1),-1);  
% B1 =  eye(N*m+1,N*m+1) + (-G_tilde)*diag(ones(N*m,1),-1);

%Initial coarse Parareal iteration at the coarse level
rhs0_coarse = [u0;zeros(N,1)];
U0_coarse = B1\rhs0_coarse;

U_fine=spline(xtrue(pick),U0_coarse,xtrue)';

% U_fine=zeros(N_fine_nodes,1);
% U_fine(1:m:end,1) =  U0_coarse;
% for i=1:m:N_fine_nodes-m
%     U_fine(i+1) = phi*U_fine(i);        
% end


% first 2-level dd preconditioning iterartion

% rhs = [u0;  phi*U(1);-G_tilde*U(1);phi*U(3);-G_tilde*U(3);phi*U(5);-G_tilde*U(5);phi*U(7);-G_tilde*U(7);phi*U(9);-G_tilde*U(9)];
%Initial fine solution at the fine level
% rhs=zeros(N_fine_nodes,1);
% rhs(1) = u0;
% for i=2:m:N_fine_nodes
%     % i = i + 2
%     rhs(i+1) = -G_tilde*U_fine(i-1) + phi/2*U_fine(i);
%     if i == N_fine_nodes - 1
%         rhs(i+1) = -G_tilde*U_fine(i-1);
%     end
% end
% U_fine = M\rhs;

% rhs=[u0];
% for i=1:2:2*N
%     rhs=[rhs; phi*U(i);-G_tilde*U(i)];
% end
% U = M\rhs;


 figure
 hold on
 plot(xtrue(pick),uF_true(pick),'-^k',xtrue(pick),U_fine(pick),'-*m');
 legend('True solution','Parareal solution at k = 0')
 
U0n = U_fine;
proNj = 1;

L2NormError = [ sqrt(sum(Dt*(uF_true(pick)-U_fine(pick)).^2))];
LInfNormError = [max(abs(uF_true(pick)-U_fine(pick)))];
SupperlinearErrorBound  = LInfNormError;
LinearErrorBound  = LInfNormError;

%%
% rhs
f = zeros(N_fine_nodes,1);
f(1) = u0;


% iteration loop
for k = 1: K-1
        % k = k + 1
        % i = i + 2
%         rhs=[u0];
%         for i=2:2:N_fine_nodes
%             rhs(i+1,1) = -G_tilde*U_fine(i-1) + phi/2*U_fine(i);
%             if i== N_fine_nodes - 1
%                 rhs(i+1) = -G_tilde*U_fine(i-1);
%             end
%         end

%         U_fine = M\rhs;
%         U_fine = U_fine + M^(-1)*(f-A*U_fine);
          U_fine = M\(M*U_fine + (f-A*U_fine));
%         U_fine = gmres(M,M*U_fine + (f-A*U_fine),10,1e-6);

       plot(xtrue,U_fine,'Color',cmap(k,:),'Marker','o');
       
%        
%        norm(uF_true(pick)-U,'inf')
       proNj(k+1) = proNj(k)*(N-k);
       L2NormError(k+1) = sqrt(sum(Dt*(uF_true(pick)-U_fine(pick)).^2));
       LInfNormError(k+1) = max(abs(uF_true(pick)-U_fine(pick)));
       SupperlinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt)))^(k)/factorial(k)*proNj(k+1)*max(abs(uF_true(pick)-U0n(pick)));
       LinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt))/(1- abs(Rz(a*Dt))))^(k)*max(abs(uF_true(pick)-U0n(pick)));
end  % iteration loop

uF_true_coarse=uF_true
U_fine

L2NormError=L2NormError'
LInfNormError=LInfNormError'
SupperlinearErrorBound=SupperlinearErrorBound'
LinearErrorBound=LinearErrorBound'

figure
semilogy(0:K-1,LInfNormError,'b--^',0:K-1,L2NormError,'g--^',0:K-1,SupperlinearErrorBound,'r*-',0:K-1,LinearErrorBound,'mx-')
legend('L^{oo}NormError','L^2NormError','Superlinear bound','Linear bound')
xlabel('k')
ylim([1e-18 1000])




















