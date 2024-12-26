close all
clear
clc


%% ---------------------------------------------------------------
% IVP: u_t = a*u_xx - b*u_x + c*u + f               in (0,L) x (0,T),
%      u(x,0) = u0(x)                    in (0,L),
%      u(0,t) = u(L,t) = 0                 t in (0,T)


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



%% -----------------------------------------------------------------


global xC xF a b c

T = 100; %  Intervall (0,T)
L = 1; % omega=(0,L)

type = 'CD'; % Heat 'H',  Convection-diffusion 'CD'


nt_interval = 100; % number of time intervals

nx_coarse = 10; % number of coarse spamtial intervals
nt_coarse = 100;% 300 %20    % number of coarse time intervals


nx_fine = 	10; % number of fine spatial intervals
nt_fine = 200;% 3000 %200 . % number of fine time intervals (nt_coarsex2 )


if strcmp(type,'H') == 1
%--------------------------------------
% % % solution x(L-x)^2*exp(-2*t);

% % %Heat
a = 3;%3
b = 0;%0.005
c = 0;%1
% 
% % %Convection-Diffusion
% % a = 3;%3
% % b = 0.005;%0.005
% % c = 1;%1
% 

%%%%--------------------------------------
% figure
% spy(P)
% figure
% spy(R'*M_tilde*R + Id - R'*R)
% figure
% spy(M_mul)

else
%--------------------------------------
% solution sin(2*pi*x)*exp(-2*t);

% % %Heat
% a = 1;      %0.5;
% b = 0;      %0.0025;
% c = 0;      %0;

% % %Convection-Diffusion
a = 0.5;
b = 0.0025;
c = 1;
% 
% nt_interval = 20; % number of time intervals
% 
% nx_coarse = 10;
% nt_coarse = 20;% 300 %20    % number of coarse time intervals
% 
% 
% nx_fine = 	10;
% nt_fine = 400;% 3000 %200 . % number of fine time intervals

end

%--------------------------------------
% % % % solution x*(L-x)*exp(-2*t);
%
% a = 3;
% b = 0.005;
% c = 1;
%
% nt_interval = 10; % number of time intervals
%
% nx_coarse = 10; % number of coarse spatial intervals
% nt_coarse = 10;% 300 %20    % number of coarse time intervals
%
%
% nx_fine = 	10; % number of fine spatial intervals
% nt_fine = 160;% 3000 %200 . % number of fine time intervals (nt_coarsex2 )


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

K = 20 % number of Parareal iterations

% 2-level domain decomposition initial
N_fine_nodes = m*nt_interval+1 % 11
% N_subdomains = m + 1% nt_interval + 1 % 6
% Id = eye(N_fine_nodes);


% plot the time domain
figure
plot(0:N_fine_nodes-1,1,'.b',0:m:N_fine_nodes-1,1,'ro')
xlim([0 N_fine_nodes-1])
ylim([0.5 1.5])

%% -----------------------------------------------------------------


xt = linspace(0,T,nt_interval+1);

cmap = hsv(K); % color map

% For the plots

% coarse points
xC = linspace(0,L,nx_coarse+1);
tC = linspace(0,T,nt_interval+1);

% fine points
xF = linspace(0,L,nx_fine+1);
tF = linspace(0,T,nt_fine+1);

cmap = hsv(K);


if strcmp(type,'H') == 1
%-------------------------------------------------------
% % exact solution 1
u_ex=@(x,t) x*(L-x)^2*exp(-2*t);%sin(2*pi*x)*exp(-2*t);%

ux0_function = @(x) (x.*(L-x).^2);%sin(2*pi*x);%

f = @(x,t) b.*(exp(-2.*t).*(L - x).^2 - x.*exp(-2.*t).*(2.*L - 2.*x)) - a.*(2.*x.*exp(-2.*t) - 2.*exp(-2.*t).*(2.*L - 2.*x)) - 2.*x.*exp(-2.*t).*(L - x).^2 - c.*x.*exp(-2.*t).*(L - x).^2;
% f = @(x,t) x.^4.*(1-x)+t.^2;

else
%-------------------------------------------------------
% % exact solution 2
u_ex=@(x,t) sin(2*pi*x)*exp(-2*t);%x*(L-x)^2*exp(-2*t);

%syms a b c x t L
% u = sin(2*pi*x)*exp(-2*t);%x*(L-x)^2*exp(-2*t)
% f = diff(u,t) - a*diff(diff(u,x),x) + b*diff(u,x) - c*u

ux0_function = @(x) sin(2*pi*x);%(x.*(L-x).^2);

f = @(x,t) 4*a*pi^2*exp(-2*t)*sin(2*pi*x) - c*exp(-2*t)*sin(2*pi*x) - 2*exp(-2*t)*sin(2*pi*x) + 2*b*pi*exp(-2*t)*cos(2*pi*x);
end
%-------------------------------------------------------
% % exact solution 3
% % % syms a b c x t L
% % % u = x*(L-x)*exp(-2*t)
% % % f = diff(u,t) - a*diff(diff(u,x),x) + b*diff(u,x) - c*u
% u_ex=@(x,t) x*(L-x)*exp(-2*t);%sin(2*pi*x)*exp(-2*t);%
%
% ux0_function = @(x) (x.*(L-x));%sin(2*pi*x);%
%
% f = @(x,t) 2*a*exp(-2*t) - b*(x.*exp(-2*t) - exp(-2*t).*(L - x)) - 2*x.*exp(-2*t).*(L - x) - c*x.*exp(-2*t).*(L - x);



%% Define fine and coarse propagators
% Backward Euler fine solver
% A_tilde
rF = a*dt/(dx^2);
kF = b*dt/(2*dx);


A_F = (kF-rF)*diag(ones(nx_coarse-2,1),1) + (1+2*rF-c*dt)*eye(nx_coarse-1) - (kF+rF)*diag(ones(nx_coarse-2,1),-1);
Id_F = eye(size(A_F));
A_F_Global =  zeros((nt_interval+1)*(nx_coarse-1));
for i =1:nt_interval+1
    %  i = i + 1
    A_F_Global((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i-1)*(nx_coarse-1)+1:i*(nx_coarse-1)) = Id_F;
    if i>1
        A_F_Global((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i-2)*(nx_coarse-1)+1:(i-1)*(nx_coarse-1)) = -A_F^-(m);
    end
end
A_tilde = A_F_Global;
condNumberA_F = cond(-A_F^-(m))
condNumberA_F_Global=cond(A_F_Global)
% Backward Euler coarse solver
% M_tilde
rC = a*dT/(dX^2);
kC = b*dT/(2*dX);

A_C = (kC-rC)*diag(ones(nx_coarse-2,1),1) + (1+2*rC-c*dT)*eye(nx_coarse-1) - (kC+rC)*diag(ones(nx_coarse-2,1),-1);
Id_C = eye(size(A_C));
A_C_Global =  zeros((nt_interval+1)*(nx_coarse-1));
for i =1:nt_interval+1
    %  i = i + 1
    A_C_Global((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i-1)*(nx_coarse-1)+1:i*(nx_coarse-1)) = Id_C;
    if i>1
        A_C_Global((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i-2)*(nx_coarse-1)+1:(i-1)*(nx_coarse-1)) = -A_C^-1;
    end
end
M_tilde = A_C_Global;
condNumberA_C = cond(-A_C^-(1))
condNumberA_C_Global = cond(A_C_Global)
%%


% Exact solution at time T = 1
U_exact=zeros(1,nx_fine+1);
U_exact_full = zeros((nt_interval+1)*(nx_fine-1),1);
for i=1:nx_fine+1
    U_exact(i)=u_ex(xF(i),T);
end

% initial solution u0
ux0 = ux0_function(xF(2:end-1))';

% % Define original A and I_F
A_Fo = (kF-rF)*diag(ones(nx_fine-2,1),1) + (1+2*rF-c*dt)*eye(nx_fine-1) - (kF+rF)*diag(ones(nx_fine-2,1),-1);
Id_Fo = eye(size(A_Fo));
A_Original =  zeros(m*nt_interval+1);
for i =1:m*nt_interval+1
    %  i = i + 1
    A_Original((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i-1)*(nx_coarse-1)+1:i*(nx_coarse-1)) = Id_Fo;
    if i>1
        A_Original((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i-2)*(nx_coarse-1)+1:(i-1)*(nx_coarse-1)) = -A_Fo^-1;
    end
end
condNumber_A_Original=cond(A_Original)



% Restriction matric for the right hand side
R_rhs =  zeros((nt_interval+1)*(nx_fine-1),(m*nt_interval+1)*(nx_fine-1));
seq_rhs = [zeros(nx_fine-1,(m*nt_interval+1)*(nx_fine-1))];
for i = 1:m
    seq_rhs(:,(i)*(nx_fine-1)+1:(i+1)*(nx_fine-1)) =  A_Fo^-(m-i);
end

for i = 1:nt_interval+1
    %  i = i + 1
    if i == 1
        R_rhs((i-1)*(nx_fine-1)+1:i*(nx_fine-1),(i-1)*(nx_fine-1)+1:i*(nx_fine-1)) = Id_Fo;
        i_temp = 0;
        
    else
        R_rhs((i-1 )*(nx_fine-1)+1:(i )*(nx_fine-1),:) = seq_rhs;
        seq_rhs = [zeros(nx_fine-1,(m)*(nx_fine-1)) seq_rhs];
        seq_rhs(:,end-(m)*(nx_fine-1)+1:end)=[];
        i_temp = i_temp + nx_fine-1;
    end
end


% source force F
F = zeros((m*nt_interval+1)*(nx_fine-1),1);
F(1:nx_fine-1) = ux0;
% i_temp = 1;
% t0 = tF(1); i = i+1
for i = 2:m*nt_interval+1
    F((i-1)*(nx_fine-1) + 1 :i*(nx_fine-1) ,:) = A_Fo^-1*dt*f(xF(2:nx_fine),tF(i))';
end

uF_seq =  FineSolver(A_F,m,nt_interval,nx_fine,F);

uF_seq_mat=[];
for i = 1:length(uF_seq)/(nx_coarse-1)
    uF_seq_mat(:,i) = [uF_seq((i-1)*(nx_coarse-1)+1:(i)*(nx_coarse-1))];
end


%     Uk_0 = uF_true;
% end
figure
hold on
plot(xF,[0;uF_seq(end-nx_fine+2:end);0])

plot(xF,U_exact)
legend('U_{FineSequential}','U_{ex}')
title('Sequentially Fine Solution and the Exact Solution')

% uF_seq = U_exact;
% return
%% ------------------------------------------------------------------
% initial coarse integration

% initial solution u0
ux0C = ux0_function(xC(2:end-1))';

% Coarse sequential propagator

Uk_0C = zeros((nt_interval+1)*(nx_coarse-1),1);
Uk_0C(1:nx_coarse-1) = ux0C;
% rhs0 = [ux0; zeros(((nt_interval+1)*(nx_coarse-1)-nx_coarse+1),1)];
%  rhs = [ux0; (A_tilde-M_tilde)*Uk_0];
% source force F
F_C = zeros((nt_interval+1)*(nx_coarse-1),1);
F_C(1:nx_coarse-1) = ux0C;
% i_temp = 1;
for i = 2:nt_interval+1
    F_C((i-1)*(nx_coarse-1) + 1 :i*(nx_coarse-1) ,:) = A_C^-1*dT*f(xC(2:nx_coarse),tC(i))';
end


% coarse space initial
indice_R_coarse  = 0:m:m*nt_interval;
R = zeros((nx_coarse-1)*length(indice_R_coarse),(m*nt_interval+1)*(nx_coarse-1));
for i = 1:length(indice_R_coarse)
    R((i-1)*(nx_coarse-1)+1:(i)*(nx_coarse-1),indice_R_coarse(i)*(nx_coarse-1)+1:indice_R_coarse(i)*(nx_coarse-1)+(nx_coarse-1)) = Id_C;
end

% coarse correction
R'*M_tilde*R;
% return
% first subdomain  Omega_1 = {0}
R1 = zeros((nx_coarse-1),(m*nt_interval+1)*(nx_coarse-1));
R1(1:nx_coarse-1,1:nx_coarse-1) = Id_C;
A1 = R1*A_Original*R1';



% define nodes in other subdomains
N_nodes_subdomain = m;

% Number of subdomains except the first subdomain
N_subdomain = m/N_nodes_subdomain*nt_interval%nt_fine/N_nodes_subdomain

indice_R_subdomains = zeros(N_subdomain,N_nodes_subdomain);
for i =1:size(indice_R_subdomains,1)
    indice_R_subdomains(i,:) =  (i-1)*N_nodes_subdomain+1:(i-1)*N_nodes_subdomain+N_nodes_subdomain;
end

% Omega_2 = {1,2} Omega_3 = {3,4}
R2_end = zeros(N_subdomain*N_nodes_subdomain*(nx_coarse-1),(N_subdomain*N_nodes_subdomain+1)*(nx_coarse-1)) ;
for i = 1:N_subdomain*N_nodes_subdomain
    R2_end((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i)*(nx_coarse-1)+1:(i+1)*(nx_coarse-1)) = Id_F;
end

R_sub = [];       % restriction matrices
A_sub = [];       % sub-matrices
for i = 2: N_subdomain+1
    R_sub{i} = R2_end(N_nodes_subdomain*(nx_coarse-1)*(i-2)+1:N_nodes_subdomain*(nx_coarse-1)*(i-2)+N_nodes_subdomain*(nx_coarse-1),:);
    A_sub{i} = R_sub{i}*A_Original*R_sub{i}';
end

% % define 2-level domain decomposition multiplicative preconditioner
P = R1'*A1*R1;
P1 = R1'*A1^-1*R1;
for i = 2:N_subdomain+1
    P = P + R_sub{i}'*A_sub{i}*R_sub{i};
    P1 = P1 + R_sub{i}'*A_sub{i}^-1*R_sub{i};
end
DDD = P^-1 - P1;
norm(DDD)
% return

%%Multiplicative
Id = eye(size(A_Original));
M_mul = P*(R'*M_tilde*R +  Id - R'*R);

%%Additive
M_add = P + (R'*M_tilde*R - R'*R);

%%% Coarse space correction
B = R'*M_tilde*R +  Id - R'*R;

% initial coarse solution
uC = M_tilde\( F_C);
uC_mat = [];


for i = 1:length(uC)/(nx_coarse-1)
    uC_mat(:,i) = [uC((i-1)*(nx_coarse-1)+1:(i)*(nx_coarse-1))];
end
% uC_mat=uC_mat';
ttrue = linspace(0,T,nt_interval*m+1);
pick = 1:m:nt_interval*m+1;
U_fine_interpolation = spline(ttrue(pick),uC_mat,ttrue);
% U_fine_interpolation = interp1(ttrue(pick),uC_mat',ttrue,'spline');

% return
U_fine = [];
for i = 1:size(U_fine_interpolation,2)
    U_fine((i-1)*(nx_coarse-1)+1:(i)*(nx_coarse-1),:) = U_fine_interpolation(:,i);
end


U_fine_mat_reduce = U_fine_interpolation(:,pick);
for i = 1: nt_interval+1
    U_fine_pick((i-1)*(nx_coarse-1)+1:(i)*(nx_coarse-1),:) =  U_fine_mat_reduce(:,i);
end






% proNj = 1;

% normInf2executiveParearealIter = [];
% normInfDefect = [norm(U_fine_pick - uF_seq,inf)];%[norm(uC(end-nx_coarse+2:end,:) - uF_seq(end-nx_coarse+2:end,:),inf)];%
%
% L2NormError = [ norm(U_fine_pick - uF_seq,2)];%[ sqrt(sum(sum(dT*(uF_seq(end-nx_coarse+2:end,:)-uC(end-nx_coarse+2:end,:)).^2)))];%
% LInfNormError = [norm(U_fine_pick - uF_seq,inf)];%[max(abs(uF_seq(end-nx_coarse+2:end,:)-uC(end-nx_coarse+2:end,:)))];%


normInf2executiveParearealIter = [];
normInfDefect = [norm(U_fine_mat_reduce' - uF_seq_mat',inf)];%[norm(uC(end-nx_coarse+2:end,:) - uF_seq(end-nx_coarse+2:end,:),inf)];%



% L2NormError_FC = [ norm(U_fine_mat_reduce' - uF_seq_mat',2)];%[ sqrt(sum(sum(dT*(uF_seq(end-nx_coarse+2:end,:)-uC(end-nx_coarse+2:end,:)).^2)))];%
L2NormError_FC = [ sqrt(sum(dX*sum(dT*(U_fine_mat_reduce' - uF_seq_mat').^2)))];%[ sqrt(sum(sum(dT*(uF_seq(end-nx_coarse+2:end,:)-uC(end-nx_coarse+2:end,:)).^2)))];%
LInfNormError_FC = [norm(U_fine_mat_reduce' - uF_seq_mat',inf)];%[max(abs(uF_seq(end-nx_coarse+2:end,:)-uC(end-nx_coarse+2:end,:)))];%

L2NormError_FCF = [ sqrt(sum(dX*sum(dT*(U_fine_mat_reduce' - uF_seq_mat').^2)))];%[ sqrt(sum(sum(dT*(uF_seq(end-nx_coarse+2:end,:)-uC(end-nx_coarse+2:end,:)).^2)))];%
LInfNormError_FCF = [norm(U_fine_mat_reduce' - uF_seq_mat',inf)];%[max(abs(uF_seq(end-nx_coarse+2:end,:)-uC(end-nx_coarse+2:end,:)))];%

L2NormError_FCFF = [ sqrt(sum(dX*sum(dT*(U_fine_mat_reduce' - uF_seq_mat').^2)))];%[ sqrt(sum(sum(dT*(uF_seq(end-nx_coarse+2:end,:)-uC(end-nx_coarse+2:end,:)).^2)))];%
LInfNormError_FCFF = [norm(U_fine_mat_reduce' - uF_seq_mat',inf)];%[max(abs(uF_seq(end-nx_coarse+2:end,:)-uC(end-nx_coarse+2:end,:)))];%

L2NormError_FCFC = [ sqrt(sum(dX*sum(dT*(U_fine_mat_reduce' - uF_seq_mat').^2)))];%[ sqrt(sum(sum(dT*(uF_seq(end-nx_coarse+2:end,:)-uC(end-nx_coarse+2:end,:)).^2)))];%
LInfNormError_FCFC = [norm(U_fine_mat_reduce' - uF_seq_mat',inf)];%[max(abs(uF_seq(end-nx_coarse+2:end,:)-uC(end-nx_coarse+2:end,:)))];%

L2NormError_FCFCF = [ sqrt(sum(dX*sum(dT*(U_fine_mat_reduce' - uF_seq_mat').^2)))];%[ sqrt(sum(sum(dT*(uF_seq(end-nx_coarse+2:end,:)-uC(end-nx_coarse+2:end,:)).^2)))];%
LInfNormError_FCFCF = [norm(U_fine_mat_reduce' - uF_seq_mat',inf)];%[max(abs(uF_seq(end-nx_coarse+2:end,:)-uC(end-nx_coarse+2:end,:)))];%

FC_boundL2 = L2NormError_FC;
FCF_boundL2 = L2NormError_FC;
FCFF_boundL2 = L2NormError_FC;
FCFC_boundL2 = L2NormError_FC;
FCFCF_boundL2 = L2NormError_FC;

%%
U0n = U_fine_mat_reduce;
% U=U_fine;


phi = A_F^-1;
phiDT = A_C^-1;
lambda = eig(phi);
mu = eig(phiDT);
U_FC = U_fine;
U_FCF = U_fine;
U_FCFF = U_fine;
U_FCFC = U_fine;
U_FCFCF = U_fine;

% iteration loop
U_fine_temp=[];
for k = 1: K
    F - A_Original*U_FC    % FC - relaxation
    [U_FC,L2NormError_FC_k,LInfNormError_FC_k] = TwolevelDDsolveConvDiff(F,A_Original,U_FC,N_subdomain,dX,dT,R1,A1,R_sub,A_sub,R,M_tilde,pick,nx_fine,uF_seq_mat,'FC')
    L2NormError_FC(k+1) = L2NormError_FC_k;
    LInfNormError_FC(k+1) = LInfNormError_FC_k;
    
    % FCF - relaxation
    [U_FCF,L2NormError_FCF_k,LInfNormError_FCF_k] = TwolevelDDsolveConvDiff(F,A_Original,U_FCF,N_subdomain,dX,dT,R1,A1,R_sub,A_sub,R,M_tilde,pick,nx_fine,uF_seq_mat,'FCF')
    L2NormError_FCF(k+1) = L2NormError_FCF_k;
    LInfNormError_FCF(k+1) = LInfNormError_FCF_k;
    
    % FCFF - relaxation
    [U_FCFF,L2NormError_FCFF_k,LInfNormError_FCFF_k] = TwolevelDDsolveConvDiff(F,A_Original,U_FCFF,N_subdomain,dX,dT,R1,A1,R_sub,A_sub,R,M_tilde,pick,nx_fine,uF_seq_mat,'FCFF')
    L2NormError_FCFF(k+1) = L2NormError_FCFF_k;
    LInfNormError_FCFF(k+1) = LInfNormError_FCFF_k;
    
    
    % FCFC - relaxation
    [U_FCFC,L2NormError_FCFC_k,LInfNormError_FCFC_k] = TwolevelDDsolveConvDiff(F,A_Original,U_FCFC,N_subdomain,dX,dT,R1,A1,R_sub,A_sub,R,M_tilde,pick,nx_fine,uF_seq_mat,'FCFC')
    L2NormError_FCFC(k+1) = L2NormError_FCFC_k;
    LInfNormError_FCFC(k+1) = LInfNormError_FCFC_k;
    
    % FCFCF - relaxation
    [U_FCFCF,L2NormError_FCFCF_k,LInfNormError_FCFCF_k] = TwolevelDDsolveConvDiff(F,A_Original,U_FCFCF,N_subdomain,dX,dT,R1,A1,R_sub,A_sub,R,M_tilde,pick,nx_fine,uF_seq_mat,'FCFCF')
    L2NormError_FCFCF(k+1) = L2NormError_FCFCF_k;
    LInfNormError_FCFCF(k+1) = LInfNormError_FCFCF_k;
    
    
    
    %     plot(xC,[0;U_fine(end-nx_coarse+2:end);0],'Color',cmap(k,:),'Marker','o');
    
    U_fine_mat_FC = [];
    for i = 1:length(U_FC)/(nx_fine-1)
        U_fine_mat_FC(:,i) = [U_FC((i-1)*(nx_fine-1)+1:(i)*(nx_fine-1))];
    end   
    U_fine_mat_reduce_FC =  U_fine_mat_FC(:,pick);
    
      U_fine_mat_FCF = [];
    for i = 1:length(U_FCF)/(nx_fine-1)
        U_fine_mat_FCF(:,i) = [U_FCF((i-1)*(nx_fine-1)+1:(i)*(nx_fine-1))];
    end   
    U_fine_mat_reduce_FCF =  U_fine_mat_FCF(:,pick);
    
      U_fine_mat_FCFF = [];
    for i = 1:length(U_FCFF)/(nx_fine-1)
        U_fine_mat_FCFF(:,i) = [U_FCFF((i-1)*(nx_fine-1)+1:(i)*(nx_fine-1))];
    end   
    U_fine_mat_reduce_FCFF =  U_fine_mat_FCFF(:,pick);
    
      U_fine_mat_FCFC = [];
    for i = 1:length(U_fine)/(nx_fine-1)
        U_fine_mat_FCFC(:,i) = [U_FCFC((i-1)*(nx_fine-1)+1:(i)*(nx_fine-1))];
    end   
    U_fine_mat_reduce_FCFC =  U_fine_mat_FCFC(:,pick);
    
      U_fine_mat_FCFCF = [];
    for i = 1:length(U_FCFCF)/(nx_fine-1)
        U_fine_mat_FCFCF(:,i) = [U_FCFCF((i-1)*(nx_fine-1)+1:(i)*(nx_fine-1))];
    end   
    U_fine_mat_reduce_FCFCF =  U_fine_mat_FCFCF(:,pick);
    
    
     FC_boundL2(k+1) = max(abs(lambda.^m - mu).*(1 - abs(mu).^((nt_interval)/m) )./ (1 - abs(mu))).^(k)*sqrt(sum(dX*sum(dT*(U_fine_mat_reduce_FC' - uC_mat').^2)));
    FCF_boundL2(k+1) = max(abs(lambda.^m - mu).*(1 - abs(mu).^((nt_interval)/m-1) )./ (1 - abs(mu)).*abs(lambda).^(m)).^(k)*sqrt(sum(dX*sum(dT*(U_fine_mat_reduce_FCF' - uC_mat').^2)));
    FCFF_boundL2(k+1) = max(abs(lambda.^m - mu).*(1 - abs(mu).^((nt_interval)/m-2) )./ (1 - abs(mu)).*abs(lambda).^(2*m)).^(k)*sqrt(sum(dX*sum(dT*(U_fine_mat_reduce_FCFF' - uC_mat').^2)));
    FCFC_boundL2(k+1) = max((lambda.^m - mu).^2.*(1 - ((nt_interval)/m ).*abs(mu).^((nt_interval)/m-1) + ((nt_interval)/m-1).*abs(mu).^((nt_interval)/m) )./ (1 - abs(mu)).^2).^(k)*sqrt(sum(dX*sum(dT*(U_fine_mat_reduce_FCFC' - uC_mat').^2)));
    FCFCF_boundL2(k+1) = max((lambda.^m - mu).^2.*(1 - ((nt_interval)/m - 1 ).*abs(mu).^((nt_interval)/m-2) + ((nt_interval)/m-2).*abs(mu).^((nt_interval)/m-1) )./ (1 - abs(mu)).^2.*abs(lambda).^m).^(k)*sqrt(sum(dX*sum(dT*(U_fine_mat_reduce_FCFCF' - uC_mat').^2)));

    
    
    
end  % iteration loop

% uF_true_coarse=uF_true(pick)
U_fine_pick;

% U_fine_mat = [];
% for i = 1:length(U_fine)/(nx_fine-1)
%     U_fine_mat(:,i) = [U_fine((i-1)*(nx_fine-1)+1:(i)*(nx_fine-1))];
% end

% U_fine_mat_reduce = U_fine_mat(:,pick);
U_plot = [zeros(nt_interval+1,1) U_fine_mat_reduce' zeros(nt_interval+1,1)];

figure
surf(xC,tC,U_plot)


% figure
% semilogy(1:K, normInf2executiveParearealIter,'-*b',0:K,normInfDefect,'-*r')
% legend('normInf2executiveParearealIter','normInfDefect')
% xlabel('k')
% title(solver)

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


% FC_boundL2 = FC_boundL2'
% FCF_boundL2 = FCF_boundL2'
% FCFF_boundL2 = FCFF_boundL2'
% FCFC_boundL2 = FCFC_boundL2'
% FCFCF_boundL2 = FCFCF_boundL2'


% figure
% semilogy(0:K,LInfNormError,'r-^',0:K,L2NormError,'b-o')
% legend('L^{\infty}NormError','L^2NormError')
% xlabel('k')
% % xlim([0 30])
% ylim([1e-16 1e-2])
% title('Error norm between the Parareal solution and the Sequential Fine solution at the end of each time-slice')

% figure
% semilogy(0:K,L2NormError,'b--o',0:K,FC_boundL2,'m-*')
% legend('L^2NormError','FC boundL2')
% xlabel('k')
% % xlim([0 30])
% ylim([1e-16 1e-2])
% title('Error norm between the Parareal solution and the Sequential Fine solution at the end of each time-slice')

% norm((R'*M_tilde*R + Id - R'*R)^-1 - (R'*M_tilde^-1*R + Id - R'*R))

%----------------------------
% %FC bound L2
figure
semilogy(0:K,L2NormError_FC,'b--^',0:K,FC_boundL2,'ro-')
legend('L^2NormError','FC boundL2')
xlabel('k')
ylim([1e-18 1])
title(['T=',num2str(T)])
%----------------------------
% %FCF bound L2
figure
semilogy(0:K,L2NormError_FCF,'b--^',0:K,FCF_boundL2,'ro-')
legend('L^2NormError','FCF boundL2')
xlabel('k')
ylim([1e-18 1])
title(['T=',num2str(T)])
%----------------------------
%FCFF bound L2
figure
semilogy(0:K,L2NormError_FCFF,'b--^',0:K,FCFF_boundL2,'ro-')
legend('L^2NormError','FCFF boundL2')
xlabel('k')
ylim([1e-18 1])
title(['T=',num2str(T)])



%----------------------------
% %FCFC bound L2
figure
semilogy(0:K,L2NormError_FCFC,'b--^',0:K,FCFC_boundL2,'ro-')
legend('L^2NormError','FCFC boundL2')
xlabel('k')
ylim([1e-18 1])
title(['T=',num2str(T)])

% %----------------------------
% % %FCFC bound L2
figure
semilogy(0:K,L2NormError_FCFCF,'b--^',0:K,FCFCF_boundL2,'ro-')
legend('L^2NormError','FCFCF boundL2')
xlabel('k')
ylim([1e-18 1])
title(['T=',num2str(T)])

% %----------------------------
% figure
% semilogy(0:K,L2NormError_FC,'b--^',0:K,L2NormError_FCF,'m--^',0:K,L2NormError_FCFF,'c--^',0:K,L2NormError_FCFC,'g--^',0:K,L2NormError_FCFCF,'r--^','LineWidth',1.5) % northeast
% legend({'L^2- FC- Two-level DD','L^2- FCF- Two-level DD','L^2- FCFF- Two-level DD','L^2- FCFC- Two-level DD','L^2- FCFCF- Two-level DD'},'Location','northeast','FontSize',12) % southwest
% xlabel('Iteration')
% ylabel('Error')
% ylim([1e-18 1e-3]) 
% title(['T=',num2str(T)])

figure
semilogy(0:K,L2NormError_FC,'b--^',0:K,L2NormError_FCF,'m--^',0:K,L2NormError_FCFF,'c--^',0:K,L2NormError_FCFCF,'r--^','LineWidth',1.5) % northeast
legend({'SC two-level additive Schwarz','SCS two-level additive Schwarz','SCS^{2} two-level additive Schwarz','S(CS)^{2} two-level additive Schwarz'},'Location','northeast','FontSize',12) % northeast southwest
xlabel('Iteration')
ylabel('Error')
ylim([1e-18 1e-1]) 
title(['T=',num2str(T)])





% figure
% semilogy(0:K,L2NormError_FC,'b-^',0:K,L2NormError_FCF,'m-^',0:K,L2NormError_FCFF,'c-^',0:K,L2NormError_FCFC,'g-^',0:K,L2NormError_FCFCF,'r-^')
% legend('L^2NormError FC-relaxation','L^2NormError FCF-relaxation','L^2NormError FCFF-relaxation','L^2NormError FCFC-relaxation','L^2NormError FCFCF-relaxation')
% xlabel('k')
% ylim([1e-18 1])
% title(['T=',num2str(T)])










% % %FC - relaxation
% FC_relax_error_propagation = (eye(n)-M1^-1*A);
% FC_relax_error_propagation_coarse = FC_relax_error_propagation(1:3:end,1:3:end)
%
% % %FCF - relaxation
% % a = (eye(n)-P^-1*A)
% % b = (eye(n)-M1^-1*A)
% % a*b - b*a
% FCF_relax_error_propagation = (eye(n)-P^-1*A)*(eye(n)-M1^-1*A);
% FCF_relax_error_propagation_coarse = FCF_relax_error_propagation(1:3:end,1:3:end)
%
%
% % %FCFF - relaxation
% FCFF_relax_error_propagation = (eye(n)-P^-1*A)*(eye(n)-P^-1*A)*(eye(n)-M1^-1*A);
% FCFF_relax_error_propagation_coarse = FCFF_relax_error_propagation(1:3:end,1:3:end)
%
%
% % % FC-FC - relaxation
% FCFC_relax_error_propagation = (eye(n)-M1^-1*A)*(eye(n)-M1^-1*A);
% FCFC_relax_error_propagation_coarse = FCFC_relax_error_propagation(1:3:end,1:3:end)
%
% % % FC-FC-F - relaxation
% FCFCF_relax_error_propagation = (eye(n)-P^-1*A)*(eye(n)-M1^-1*A)^2;
% FCFCF_relax_error_propagation_coarse = FCFCF_relax_error_propagation(1:3:end,1:3:end)

% phi = -A_F^-1;
% phiDT = -A_C^-1;
% eig_F = eig(phi)
% eig_C = eig(phiDT)
%
% L2NormError_FC = L2NormError'
% FC_boundL2 = FC_boundL2'
% FC_bound = max(abs(eig(phi).^m - eig(phiDT)).*(1 - abs(eig(phiDT)).^((nt_interval-1)/m) )./ (1 - eig(phiDT)))
%
% FCF_bound = max(abs(eig(phi).^m - eig(phiDT)).*(1 - abs(eig(phiDT)).^((nt_interval-1)/m) )./ (1 - eig(phiDT)).*abs(eig(phi)).^m)






% test commute
a1 = 2*(phi^m-phiDT)*phiDT*(phi^m-phiDT);
a2 = 2*phiDT*(phi^m-phiDT)^2;









