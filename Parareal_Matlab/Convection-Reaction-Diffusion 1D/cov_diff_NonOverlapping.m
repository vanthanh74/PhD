close all
clear
clc


%% ---------------------------------------------------------------
% IVP: u_t = a*u_xx - b*u_x + c*u + f               in (0,L) x (0,T),
%      u(x,0) = u0(x)                    in (0,L),
%      u(0,t) = u(L,t) = 0                 t in (0,T)  

% non-overlapping subdomains 2-level domain decomposition multiplicative preconditioner

%% -----------------------------------------------------------------
% INPUT PARAMETERS:
%   * nt_coarse: number of coarse intervals in time
%   * nx_coarse: number of coarse intervals in space
%   * nt_fine: number of fine intervals in time
%   * nx_fine: number of fine intervals  in space
%   * m: number of fine time steps on each coarse time step 
%   * K: number of Parareal iterations (2 <= K <= n_coarse+1)
%        (--> K sequential coarse sweeps and K-1 parallel fine sweeps)
% OUTPUT:
%   * LInfinityErrorNorm

global xC xF a b c

T = 1; %  Intervall (0,T)
L = 1; % omega=(0,L)

a = 3;
b = 0.0005;
c = 1;

nt_interval = 2; % number of time intervals

nx_coarse = 5;
nt_coarse = 2;% 300 %20    % number of coarse time intervals


nx_fine = 	5;
nt_fine = 4;% 3000 %200 . % number of fine time intervals


% solver ='BackwardEuler';
solver ='BackwardEuler';

dx = L/nx_fine;   % fine spatial discretization steps
if strcmp(solver,'BackwardEuler')==1
    dt = T/nt_fine;   % fine temporal discretization steps . % dx^2/2 for Runge-Kutta
elseif strcmp(solver,'Runge-Kutta4')==1
    dt = dx^2/6;
end

DT = T/nt_interval;  % reference temporal discretization steps

dX = L/nx_coarse; % coarse spatial discretization steps
dT = T/nt_coarse; % coarse temporal discretization steps

M = DT/dT  % number of coarse time steps on each time interval 
m = DT/dt  % number of fine time steps on each time interval 

r = dT/dX^2;
r1 = dt/dx^2;

K = 20    % number of Parareal iterations

%% -----------------------------------------------------------------

% global xC xF a b c

% stability function (  y^n+1 = R(z)y^n  )
% Rz = @(z) 1/(1-z);     % Backward Euler
% Rz = @(z) 1+z;       % Forward Euler

%% -----------------------------------------------------------------


% % initial condition
% u0 = 1;
% 
% % coarse propagator
% N = 2;     %   number of coarse time intervals
% T = 1;      %   T = N*Dt  
% Dt = T/N;   %   coarse timesteps

xt = linspace(0,T,nt_interval+1);

% fine propagator
% m = 2;    % number of fine time steps in each coarse time interval
% dt = Dt/m;
% 
% K = 7 % number of Parareal iterations

cmap = hsv(K); % color map

% 2-level domain decomposition initial
N_fine_nodes = m*nt_interval+1; % 11
N_subdomains = m+1; % 6
Id = eye(N_fine_nodes);

% plot the time domain
figure 
plot(0:N_fine_nodes-1,1,'.b',0:m:N_fine_nodes-1,1,'ro')
xlim([0 N_fine_nodes-1])
ylim([0.5 1.5])


% For the plots


% coarse points
xC = linspace(0,L,nx_coarse+1);
tC = linspace(0,T,nt_interval+1);

% fine points
xF = linspace(0,L,nx_fine+1);
tF = linspace(0,T,nt_fine+1);

cmap = hsv(K);

% exact solution
u_ex=@(x,t) x*(L-x)^2*exp(-2*t);

ux0_function = @(x) (x.*(L-x).^2);

f = @(x,t) b*(exp(-2*t)*(L - x)^2 - x*exp(-2*t)*(2*L - 2*x)) - a*(2*x*exp(-2*t) - 2*exp(-2*t)*(2*L - 2*x)) - 2*x*exp(-2*t)*(L - x)^2 - c*x*exp(-2*t)*(L - x)^2;

% initial solution u0
ux0 = ux0_function(xF);

%% Define find and coarse propagators
% Backward Euler fine solver
rF = a*dt/(dx^2);
kF = b*dt/(2*dx);

A_F = (kF-rF)*diag(ones(nx_fine-2,1),1) + (1+2*rF-c*dt)*eye(nx_fine-1) - (kF+rF)*diag(ones(nx_fine-2,1),-1);
Id_F = eye(size(A_F));
A_F_Global =  zeros(nt_interval*(nx_fine-1));
for i =1:N_fine_nodes
    %  i = i + 1
    A_F_Global((i-1)*(nx_fine-1)+1:i*(nx_fine-1),(i-1)*(nx_fine-1)+1:i*(nx_fine-1)) = Id_F;
    if i>1
        A_F_Global((i-1)*(nx_fine-1)+1:i*(nx_fine-1),(i-2)*(nx_fine-1)+1:(i-1)*(nx_fine-1)) = -A_F^-1; 
    end
end

% Backward Euler coarse solver
rC = a*dT/(dX^2);
kC = b*dT/(2*dX);

A_C = (kC-rC)*diag(ones(nx_coarse-2,1),1) + (1+2*rC-c*dT)*eye(nx_coarse-1) - (kC+rC)*diag(ones(nx_coarse-2,1),-1);
Id_C = eye(size(A_C));
A_C_Global =  zeros(nt_interval*(nx_coarse-1));
for i =1:nx_coarse
    %  i = i + 1
    A_C_Global((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i-1)*(nx_coarse-1)+1:i*(nx_coarse-1)) = Id_C;
    if i>1
        A_C_Global((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i-2)*(nx_coarse-1)+1:(i-1)*(nx_coarse-1)) = -A_C^-1; 
    end
end



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
R2_end = zeros((N_subdomains-2)*N_nodes_subdomain,N_fine_nodes) ;
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
M = P*(R'*B*R +  Id-R'*R)


% preconditioned stationary iteration
% U^(k+1) = U^(k) + M^(-1)*( f-A*U^(k)) )




%% ------------------------------------------------------------------

% initial coarse integration

B1 =  eye(N+1,N+1) + (-G_tilde)*diag(ones(N,1),-1);  
% B1 =  eye(N*m+1,N*m+1) + (-G_tilde)*diag(ones(N*m,1),-1);

%Initial coarse Parareal iteration at the coarse level
rhs0_coarse = [u0;zeros(N,1)];
U0_coarse = B1\rhs0_coarse;

U_fine=spline(xtrue(pick),U0_coarse,xtrue)';

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
 plot(xtrue,uF_true,'-^k',xtrue,U_fine,'-*m');
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
%         rhs=[u0];
%         for i=1:2:2*N
%             rhs=[rhs; phi*U_fine(i);-G_tilde*U_fine(i)];
%         end
%         U_fine = M\rhs;
%         U_fine = U_fine + M^(-1)*(f-A*U_fine);
        U_fine = M\(M*U_fine + (f-A*U_fine));
%         U_fine = gmres(M,M*U_fine + (f-A*U_fine),10,1e-6);
       plot(xtrue(pick),U_fine(pick),'Color',cmap(k,:),'Marker','o');
       
%        
%        norm(uF_true(pick)-U,'inf')
       proNj(k+1) = proNj(k)*(N-k);
       L2NormError(k+1) = sqrt(sum(Dt*(uF_true(pick)-U_fine(pick)).^2));
       LInfNormError(k+1) = max(abs(uF_true(pick)-U_fine(pick)));
       SupperlinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt)))^(k)/factorial(k)*proNj(k+1)*max(abs(uF_true(pick)-U0n(pick)));
       LinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt))/(1- abs(Rz(a*Dt))))^(k)*max(abs(uF_true(pick)-U0n(pick)));
end  % iteration loop

uF_true=uF_true
U_fine

L2NormError=L2NormError'
LInfNormError=LInfNormError'
SupperlinearErrorBound=SupperlinearErrorBound'
LinearErrorBound=LinearErrorBound'

figure
semilogy(0:K-1,LInfNormError,'b--^',0:K-1,L2NormError,'g--^',0:K-1,SupperlinearErrorBound,'r*-',0:K-1,LinearErrorBound,'mx-')
legend('L^{oo}NormError','L^2NormError','Superlinear bound','Linear bound')
xlabel('k')
ylim([1e-18 1])




















