function [x,flag,relres,iter,resvec,L2NormErrorGMRES,LInfNormErrorGMRES] = mygmresConvDiff(A,b,restart,tol,maxit,M1,M2,x,type,varargin)
%GMRES   Generalized Minimum Residual Method.                                      (MinvAfun,rhs,[],[],K,[],[],U0_fine,type)
%   X = GMRES(A,B) attempts to solve the system of linear equations A*X = B
%   for X.  The N-by-N coefficient matrix A must be square and the right
%   hand side column vector B must have length N. This uses the unrestarted
%   method with MIN(N,10) total iterations.
%
%   X = GMRES(AFUN,B) accepts a function handle AFUN instead of the matrix
%   A. AFUN(X) accepts a vector input X and returns the matrix-vector
%   product A*X. In all of the following syntaxes, you can replace A by
%   AFUN.
%
%   X = GMRES(A,B,RESTART) restarts the method every RESTART iterations.
%   If RESTART is N or [] then GMRES uses the unrestarted method as above.
%
%   X = GMRES(A,B,RESTART,TOL) specifies the tolerance of the method.  If
%   TOL is [] then GMRES uses the default, 1e-6.
%
%   X = GMRES(A,B,RESTART,TOL,MAXIT) specifies the maximum number of outer
%   iterations. Note: the total number of iterations is RESTART*MAXIT. If
%   MAXIT is [] then GMRES uses the default, MIN(N/RESTART,10). If RESTART
%   is N or [] then the total number of iterations is MAXIT.
%
%   X = GMRES(A,B,RESTART,TOL,MAXIT,M) and
%   X = GMRES(A,B,RESTART,TOL,MAXIT,M1,M2) use preconditioner M or M=M1*M2
%   and effectively solve the system inv(M)*A*X = inv(M)*B for X. If M is
%   [] then a preconditioner is not applied.  M may be a function handle
%   returning M\X.
%
%   X = GMRES(A,B,RESTART,TOL,MAXIT,M1,M2,X0) specifies the first initial
%   guess. If X0 is [] then GMRES uses the default, an all zero vector.
%
%   [X,FLAG] = GMRES(A,B,...) also returns a convergence FLAG:
%    0 GMRES converged to the desired tolerance TOL within MAXIT iterations.
%    1 GMRES iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 GMRES stagnated (two consecutive iterates were the same).
%
%   [X,FLAG,RELRES] = GMRES(A,B,...) also returns the relative residual
%   NORM(B-A*X)/NORM(B). If FLAG is 0, then RELRES <= TOL. Note with
%   preconditioners M1,M2, the residual is NORM(M2\(M1\(B-A*X))).
%
%   [X,FLAG,RELRES,ITER] = GMRES(A,B,...) also returns both the outer and
%   inner iteration numbers at which X was computed: 0 <= ITER(1) <= MAXIT
%   and 0 <= ITER(2) <= RESTART.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = GMRES(A,B,...) also returns a vector of
%   the residual norms at each inner iteration, including NORM(B-A*X0).
%   Note with preconditioners M1,M2, the residual is NORM(M2\(M1\(B-A*X))).
%
%   Example:
%      n = 21; A = gallery('wilk',n);  b = sum(A,2);
%      tol = 1e-12;  maxit = 15; M = diag([10:-1:1 1 1:10]);
%      x = gmres(A,b,10,tol,maxit,M);
%   Or, use this matrix-vector product function
%      %-----------------------------------------------------------------%
%      function y = afun(x,n)
%      y = [0; x(1:n-1)] + [((n-1)/2:-1:0)'; (1:(n-1)/2)'].*x+[x(2:n); 0];
%      %-----------------------------------------------------------------%
%   and this preconditioner backsolve function
%      %------------------------------------------%
%      function y = mfun(r,n)
%      y = r ./ [((n-1)/2:-1:1)'; 1; (1:(n-1)/2)'];
%      %------------------------------------------%
%   as inputs to GMRES:
%      x1 = gmres(@(x)afun(x,n),b,10,tol,maxit,@(x)mfun(x,n));
%
%   Class support for inputs A,B,M1,M2,X0 and the output of AFUN:
%      float: double
%
%   See also BICG, BICGSTAB, BICGSTABL, CGS, LSQR, MINRES, PCG, QMR, SYMMLQ,
%   TFQMR, ILU, FUNCTION_HANDLE.

%   References
%   H.F. Walker, "Implementation of the GMRES Method Using Householder
%   Transformations", SIAM J. Sci. Comp. Vol 9. No 1. January 1988.

%   Copyright 1984-2015 The MathWorks, Inc.




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


global xC xF T L a bb  c nt_interval nx_coarse nt_coarse nx_fine nt_fine 

% T = 1; %  Intervall (0,T)
% L = 1; % omega=(0,L)

%--------------------------------------
% % % solution x(L-x)^2*exp(-2*t);

% %Heat
% a = 3;%3
% bb = 0.05;%0.005
% c = 0;%1


%Convection-Diffusion
% a = 3;
% bb = 0.005;
% c = 1;

% nt_interval = 20; % number of time intervals
% 
% nx_coarse = 10; % number of coarse spamtial intervals
% nt_coarse = 20;% 300 %20    % number of coarse time intervals
% 
% 
% nx_fine = 	10; % number of fine spatial intervals
% nt_fine = 40;% 3000 %200 . % number of fine time intervals (nt_coarsex2 )

%%%%--------------------------------------
% figure
% spy(P)
% figure
% spy(R'*M_tilde*R + Id - R'*R)
% figure
% spy(M_mul)

%--------------------------------------
% solution sin(2*pi*x)*exp(-2*t);

% a = 0.5;      %0.5;
% b = 0;      %0.0025;
% c = 0;      %0;
% 
% nt_interval = 20; % number of time intervals
% 
% nx_coarse = 10;
% nt_coarse = 20;% 300 %20    % number of coarse time intervals
% 
% 
% nx_fine = 	10;
% nt_fine = 40;% 3000 %200 . % number of fine time intervals



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


% % plot the time domain
% figure
% plot(0:N_fine_nodes-1,1,'.b',0:m:N_fine_nodes-1,1,'ro')
% xlim([0 N_fine_nodes-1])
% ylim([0.5 1.5])

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

f = @(x,t) bb.*(exp(-2.*t).*(L - x).^2 - x.*exp(-2.*t).*(2.*L - 2.*x)) - a.*(2.*x.*exp(-2.*t) - 2.*exp(-2.*t).*(2.*L - 2.*x)) - 2.*x.*exp(-2.*t).*(L - x).^2 - c.*x.*exp(-2.*t).*(L - x).^2;
% % f = @(x,t) x.^4.*(1-x)+t.^2;
else

%-------------------------------------------------------
% % exact solution 2
u_ex=@(x,t) sin(2*pi*x)*exp(-2*t);%x*(L-x)^2*exp(-2*t);

%syms a b c x t L
% u = sin(2*pi*x)*exp(-2*t);%x*(L-x)^2*exp(-2*t)
% f = diff(u,t) - a*diff(diff(u,x),x) + b*diff(u,x) - c*u

ux0_function = @(x) sin(2*pi*x);%(x.*(L-x).^2);

f = @(x,t) 4*a*pi^2*exp(-2*t)*sin(2*pi*x) - c*exp(-2*t)*sin(2*pi*x) - 2*exp(-2*t)*sin(2*pi*x) + 2*bb*pi*exp(-2*t)*cos(2*pi*x);

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
kF = bb*dt/(2*dx);


A_F = (kF-rF)*diag(ones(nx_coarse-2,1),1) + (1+2*rF-c*dt)*eye(nx_coarse-1) - (kF+rF)*diag(ones(nx_coarse-2,1),-1);
Id_F = eye(size(A_F));
% A_F_Global =  zeros((nt_interval+1)*(nx_coarse-1));
% for i =1:nt_interval+1
%     %  i = i + 1
%     A_F_Global((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i-1)*(nx_coarse-1)+1:i*(nx_coarse-1)) = Id_F;
%     if i>1
%         A_F_Global((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i-2)*(nx_coarse-1)+1:(i-1)*(nx_coarse-1)) = -A_F^-(m);
%     end
% end
% A_tilde = A_F_Global;
% condNumberA_F = cond(-A_F^-(m))
% condNumberA_F_Global=cond(A_F_Global)
% Backward Euler coarse solver
% M_tilde
rC = a*dT/(dX^2);
kC = bb*dT/(2*dX);

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
% condNumberA_C = cond(-A_C^-(1))
% condNumberA_C_Global = cond(A_C_Global)
%%
% Reduced system
% F_tilde = (1-a*dt)^-m;  % fine propagator
% G_tilde = (1-a*Dt)^-1;  % coarse propagator
%
% Uk_0 = zeros(N,1);
% % Uk_0(1) = u0;
% rhs0 = [u0; zeros(m*N,1)];
% rhs = [u0; (F_tilde-G_tilde)*Uk_0];
%
% A = eye(m*N+1,m*N+1) + (-(1-a*dt)^-1)*diag(ones(m*N,1),-1);  % original problem matrix
%
% A_tilde = eye(N+1,N+1) + (-F_tilde)*diag(ones(N,1),-1);      % reduced matrix
%
% M_tilde = eye(N+1,N+1) + (-G_tilde)*diag(ones(N,1),-1);      % approximate reduced matrix

% fine solution

%  xtrue = linspace(0,T,nx_coarse*m+1);
%  ttrue = linspace()
%  pick = 1:m:nt_interval*m+1;

% Exact solution at time T = 1
U_exact=zeros(1,nx_fine+1);
U_exact_full = zeros((nt_interval+1)*(nx_fine-1),1);
for i=1:nx_fine+1
    U_exact(i)=u_ex(xF(i),T);
end

% initial solution u0
ux0 = ux0_function(xF(2:end-1))';

% Fine sequential propagator
% Uk_0 = zeros((nt_interval+1)*(nx_fine-1),1);
% Uk_0(1:nx_fine-1) = ux0;
% rhs0 = [ux0; zeros(((nt_interval+1)*(nx_fine-1)-nx_fine+1),1)];
%  rhs = [ux0; (A_tilde-M_tilde)*Uk_0];

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
% condNumber_A_Original=cond(A_Original)



% % Restriction matric for the right hand side
% R_rhs =  zeros((nt_interval+1)*(nx_fine-1),(m*nt_interval+1)*(nx_fine-1));
% seq_rhs = [zeros(nx_fine-1,(m*nt_interval+1)*(nx_fine-1))];
% for i = 1:m
%     seq_rhs(:,(i)*(nx_fine-1)+1:(i+1)*(nx_fine-1)) =  A_Fo^-(m-i);
% end
% 
% for i = 1:nt_interval+1
%     %  i = i + 1
%     if i == 1
%         R_rhs((i-1)*(nx_fine-1)+1:i*(nx_fine-1),(i-1)*(nx_fine-1)+1:i*(nx_fine-1)) = Id_Fo;
%         i_temp = 0;
%         
%     else
%         %         R_rhs((i-1)*(nx_fine-1)+1:i*(nx_fine-1),(i-1+i_temp)*(nx_fine-1)+1:(i+i_temp)*(nx_fine-1)) = A_F^-1;
%         %         R_rhs((i-1)*(nx_fine-1)+1:i*(nx_fine-1),(i+i_temp)*(nx_fine-1)+1:(i+1+i_temp)*(nx_fine-1)) = Id_F;
%         %         i_temp = i_temp + 1;
%         
%         R_rhs((i-1 )*(nx_fine-1)+1:(i )*(nx_fine-1),:) = seq_rhs;
%         seq_rhs = [zeros(nx_fine-1,(m)*(nx_fine-1)) seq_rhs];
%         seq_rhs(:,end-(m)*(nx_fine-1)+1:end)=[];
%         i_temp = i_temp + nx_fine-1;
%     end
% end


% source force F
F = zeros((m*nt_interval+1)*(nx_fine-1),1);
F(1:nx_fine-1) = ux0;
% i_temp = 1;
% t0 = tF(1); i = i+1
for i = 2:m*nt_interval+1
    F((i-1)*(nx_fine-1) + 1 :i*(nx_fine-1) ,:) = A_Fo^-1*dt*f(xF(2:nx_fine),tF(i))';
end

% uF_seq = A_tilde\( R_rhs*F);
% Right_hand_side = R_rhs*F;
% u0 = Right_hand_side(1:nx_fine-1,1);
uF_seq =  FineSolver(A_F,m,nt_interval,nx_fine,F);

uF_seq_mat=[];
for i = 1:length(uF_seq)/(nx_coarse-1)
    uF_seq_mat(:,i) = [uF_seq((i-1)*(nx_coarse-1)+1:(i)*(nx_coarse-1))];
end


%     Uk_0 = uF_true;
% end
% figure
% hold on
% plot(xF,[0;uF_seq(end-nx_fine+2:end);0])
% 
% plot(xF,U_exact)
% legend('U_{FineSequential}','U_{ex}')
% title('Sequentially Fine Solution and the Exact Solution')

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

% figure
% hold on
% uC = M_tilde\( F_C);
% plot(xC,[0;uC(end-nx_coarse+2:end);0])

% plot(xF,U_exact)
% legend('U_{0}','U_{ex}')
% title('Parareal Solution and the Exact Solution')

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



% A_Fo = (kF-rF)*diag(ones(nx_fine-2,1),1) + (1+2*rF-c*dt)*eye(nx_fine-1) - (kF+rF)*diag(ones(nx_fine-2,1),-1);
% Id_Fo = eye(size(A_Fo));
% A_Original =  zeros(m*(nt_interval+1)*(nx_coarse-1));



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


% P_sum_inverse = R1'*A1^-1*R1;
% for i = 2:nt_interval+1
%     P_sum_inverse = P_sum_inverse + R_sub{i}'*A_sub{i}^-1*R_sub{i};
% end



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

U0_fine = U_fine;
% % % interp1
% for i = 1:size(U_fine_interpolation,1)
%     U_fine((i-1)*(nx_coarse-1)+1:(i)*(nx_coarse-1),:) = U_fine_interpolation(i,:)';
% end

% U_fine_mat = [];
% for i = 1:length(U_fine)/(nx_fine-1)
%     U_fine_mat(:,i) = [U_fine((i-1)*(nx_fine-1)+1:(i)*(nx_fine-1))];
% end

U_fine_mat_reduce = U_fine_interpolation(:,pick);
for i = 1: nt_interval+1
    U_fine_pick((i-1)*(nx_coarse-1)+1:(i)*(nx_coarse-1),:) =  U_fine_mat_reduce(:,i);
end

% % % interp1
% U_fine_mat_reduce = U_fine_interpolation(pick,:);
% for i = 1: nt_interval+1
%     U_fine_pick((i-1)*(nx_coarse-1)+1:(i)*(nx_coarse-1),:) =  U_fine_mat_reduce(i,:)';
% end

% figure
% surf(xC(2:end-1),tC,uC_mat')
% % hold on

% figure
% surf(U_fine_mat)
% return






% proNj = 1;

% normInf2executiveParearealIter = [];
% normInfDefect = [norm(U_fine_pick - uF_seq,inf)];%[norm(uC(end-nx_coarse+2:end,:) - uF_seq(end-nx_coarse+2:end,:),inf)];%
% 
% L2NormError = [ norm(U_fine_pick - uF_seq,2)];%[ sqrt(sum(sum(dT*(uF_seq(end-nx_coarse+2:end,:)-uC(end-nx_coarse+2:end,:)).^2)))];%
% LInfNormError = [norm(U_fine_pick - uF_seq,inf)];%[max(abs(uF_seq(end-nx_coarse+2:end,:)-uC(end-nx_coarse+2:end,:)))];%


normInf2executiveParearealIter = [];
normInfDefect = [norm(U_fine_mat_reduce' - uF_seq_mat',inf)];%[norm(uC(end-nx_coarse+2:end,:) - uF_seq(end-nx_coarse+2:end,:),inf)];%

L2NormErrorGMRES = [ sqrt(sum(dX*sum(dT*(U_fine_mat_reduce' - uF_seq_mat').^2)))];%[ sqrt(sum(sum(dT*(uF_seq(end-nx_coarse+2:end,:)-uC(end-nx_coarse+2:end,:)).^2)))];%
LInfNormErrorGMRES = [norm(U_fine_mat_reduce' - uF_seq_mat',inf)];%[max(abs(uF_seq(end-nx_coarse+2:end,:)-uC(end-nx_coarse+2:end,:)))];%


% SupperlinearErrorBound  = LInfNormError;
% LinearErrorBound  = LInfNormError;

%%
U0n = U_fine_mat_reduce;
res_Parareal=[];
% iteration loop
U_fine_temp=[norm(F-A_Original*U0_fine)];































if (nargin < 2)
    error(message('MATLAB:gmres:NumInputs'));
end

% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = iterchk(A);
if strcmp(atype,'matrix')
    % Check matrix and right hand side vector inputs have appropriate sizes
    [mm,n] = size(A);
    if (mm ~= n)
        error(message('MATLAB:gmres:SquareMatrix'));
    end
    if ~isequal(size(b),[mm,1])
        error(message('MATLAB:gmres:VectorSize', mm));
    end
else
    mm = size(b,1);
    n = mm;
    if ~iscolumn(b)
        error(message('MATLAB:gmres:Vector'));
    end
end

% Assign default values to unspecified parameters
if (nargin < 3) || isempty(restart) || (restart == n)
    restarted = false;
else
    restarted = true;
end
if (nargin < 4) || isempty(tol)
    tol = 1e-16;
end
warned = 0;
if tol < eps
    warning(message('MATLAB:gmres:tooSmallTolerance'));
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning(message('MATLAB:gmres:tooBigTolerance'));
    warned = 1;
    tol = 1-eps;
end
% if (nargin < 5) || isempty(maxit)
%     if restarted
%         maxit = min(ceil(n/restart),10);
%     else
%         maxit = K;%min(n,10);
%     end
% end
maxit = K;
inner = maxit;
 outer = 1;
% if restarted
%     outer = maxit;
%     if restart > n
%         warning(message('MATLAB:gmres:tooManyInnerItsRestart',restart, n));
%         restart = n;
%     end
%     inner = restart;
% else
%     outer = 1;
%     if maxit > n
%         warning(message('MATLAB:gmres:tooManyInnerItsMaxit',maxit, n));
%         maxit = n;
%     end
%     inner = maxit;
% end

% Check for all zero right hand side vector => all zero solution
n2b = norm(b);                   % Norm of rhs vector, b
if (n2b == 0)                    % if    rhs vector is all zeros
    x = zeros(n,1);              % then  solution is all zeros
    flag = 0;                    % a valid solution has been obtained
    relres = 0;                  % the relative residual is actually 0/0
    iter = [0 0];                % no iterations need be performed
    resvec = 0;                  % resvec(1) = norm(b-A*x) = norm(0)
    if (nargout < 2)
        itermsg('gmres',tol,maxit,0,flag,iter,NaN);
    end
    return
end

if ((nargin >= 6) && ~isempty(M1))
    existM1 = 1;
    [m1type,m1fun,m1fcnstr] = iterchk(M1);
    if strcmp(m1type,'matrix')
        if ~isequal(size(M1),[mm,mm])
            error(message('MATLAB:gmres:PreConditioner1Size', mm));
        end
    end
else
    existM1 = 0;
    m1type = 'matrix';
end

if ((nargin >= 7) && ~isempty(M2))
    existM2 = 1;
    [m2type,m2fun,m2fcnstr] = iterchk(M2);
    if strcmp(m2type,'matrix')
        if ~isequal(size(M2),[mm,mm])
            error(message('MATLAB:gmres:PreConditioner2Size', mm));
        end
    end
else
    existM2 = 0;
    m2type = 'matrix';
end

if ((nargin >= 8) && ~isempty(x))
    if ~isequal(size(x),[n,1])
        error(message('MATLAB:gmres:XoSize', n));
    end
else
    x = zeros(n,1);
end

if ((nargin > 8) && strcmp(atype,'matrix') && ...
        strcmp(m1type,'matrix') && strcmp(m2type,'matrix'))
    error(message('MATLAB:gmres:TooManyInputs'));
end

% Set up for the method
flag = 1;
xmin = x;                        % Iterate which has minimal residual so far
imin = 0;                        % "Outer" iteration at which xmin was computed
jmin = 0;                        % "Inner" iteration at which xmin was computed
tolb = tol * n2b;                % Relative tolerance
evalxm = 0;
stag = 0;
moresteps = 0;
maxmsteps = min([floor(n/50),5,n-maxit]);
maxstagsteps = 3;
minupdated = 0;

x0iszero = (norm(x) == 0);
r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
normr = norm(r);                 % Norm of initial residual
if (normr <= tolb)               % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2b;
    iter = [0 0];
    resvec = normr;
    if (nargout < 2)
        itermsg('gmres',tol,maxit,[0 0],flag,iter,relres);
    end
    return
end
minv_b = b;

if existM1
    r = iterapp('mldivide',m1fun,m1type,m1fcnstr,r,varargin{:});
    if ~x0iszero
        minv_b = iterapp('mldivide',m1fun,m1type,m1fcnstr,b,varargin{:});
    else
        minv_b = r;
    end
    if ~all(isfinite(r)) || ~all(isfinite(minv_b))
        flag = 2;
        x = xmin;
        relres = normr / n2b;
        iter = [0 0];
        resvec = normr;
        return
    end
end

if existM2
    r = iterapp('mldivide',m2fun,m2type,m2fcnstr,r,varargin{:});
    if ~x0iszero
        minv_b = iterapp('mldivide',m2fun,m2type,m2fcnstr,minv_b,varargin{:});
    else
        minv_b = r;
    end
    if ~all(isfinite(r)) || ~all(isfinite(minv_b))
        flag = 2;
        x = xmin;
        relres = normr / n2b;
        iter = [0 0];
        resvec = normr;
        return
    end
end

normr = norm(r);                 % norm of the preconditioned residual
n2minv_b = norm(minv_b);         % norm of the preconditioned rhs
clear minv_b;
tolb = tol * n2minv_b;
if (normr <= tolb)               % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2minv_b;
    iter = [0 0];
    resvec = n2minv_b;
    if (nargout < 2)
        itermsg('gmres',tol,maxit,[0 0],flag,iter,relres);
    end
    return
end

resvec = zeros(inner*outer+1,1);  % Preallocate vector for norm of residuals
resvec(1) = 1;%Relative residual%norm(F-A_Original*U0_fine);  %              % resvec(1) = norm(b-A*x0)
normrmin = normr;                 % Norm of residual from xmin

%  Preallocate J to hold the Given's rotation constants.
J = zeros(2,inner);

U = zeros(n,inner);
R = zeros(inner,inner);
w = zeros(inner+1,1);

normr0=normr;
for outiter = 1 : outer
    %  Construct u for Householder reflector.
    %  u = r + sign(r(1))*||r||*e1
    u = r;
    normr = norm(r);
    beta = scalarsign(r(1))*normr;
    u(1) = u(1) + beta;
    u = u / norm(u);
    
    U(:,1) = u;
    
    %  Apply Householder projection to r.
    %  w = r - 2*u*u'*r;
    w(1) = -beta;
    
    for initer = 1 : inner
        %  Form P1*P2*P3...Pj*ej.
        %  v = Pj*ej = ej - 2*u*u'*ej
        v = -2*(u(initer)')*u;
        v(initer) = v(initer) + 1;
        %  v = P1*P2*...Pjm1*(Pj*ej)
        for k = (initer-1):-1:1
            Utemp = U(:,k);
            v = v - Utemp*(2*(Utemp'*v));
        end
        %  Explicitly normalize v to reduce the effects of round-off.
        v = v/norm(v);
        
        %  Apply A to v.
        v = iterapp('mtimes',afun,atype,afcnstr,v,varargin{:});
        %  Apply Preconditioner.
        if existM1
            v = iterapp('mldivide',m1fun,m1type,m1fcnstr,v,varargin{:});
            if ~all(isfinite(v))
                flag = 2;
                break
            end
        end
        
        if existM2
            v = iterapp('mldivide',m2fun,m2type,m2fcnstr,v,varargin{:});
            if ~all(isfinite(v))
                flag = 2;
                break
            end
        end
        %  Form Pj*Pj-1*...P1*Av.
        for k = 1:initer
            Utemp = U(:,k);
            v = v - Utemp*(2*(Utemp'*v));
        end
        
        %  Determine Pj+1.
        if (initer ~= length(v))
            %  Construct u for Householder reflector Pj+1.
            u = v;
            u(1:initer) = 0;
            alpha = norm(u);
            if (alpha ~= 0)
                alpha = scalarsign(v(initer+1))*alpha;
                %  u = v(initer+1:end) +
                %        sign(v(initer+1))*||v(initer+1:end)||*e_{initer+1)
                u(initer+1) = u(initer+1) + alpha;
                u = u / norm(u);
                U(:,initer+1) = u;
                
                %  Apply Pj+1 to v.
                %  v = v - 2*u*(u'*v);
                v(initer+2:end) = 0;
                v(initer+1) = -alpha;
            end
        end
        
        %  Apply Given's rotations to the newly formed v.
        for colJ = 1:initer-1
            tmpv = v(colJ);
            v(colJ)   = conj(J(1,colJ))*v(colJ) + conj(J(2,colJ))*v(colJ+1);
            v(colJ+1) = -J(2,colJ)*tmpv + J(1,colJ)*v(colJ+1);
        end
        
        %  Compute Given's rotation Jm.
        if ~(initer==length(v))
            rho = norm(v(initer:initer+1));
            J(:,initer) = v(initer:initer+1)./rho;
            w(initer+1) = -J(2,initer).*w(initer);
            w(initer) = conj(J(1,initer)).*w(initer);
            v(initer) = rho;
            v(initer+1) = 0;
        end
        
        R(:,initer) = v(1:inner);
        
        normr = abs(w(initer+1));
        resvec((outiter-1)*inner+initer+1) = normr/normr0;
        normr_act = normr;
        
%         if (normr <= tolb || stag >= maxstagsteps || moresteps)
            if evalxm == 0
                ytmp = R(1:initer,1:initer) \ w(1:initer);
                additive = U(:,initer)*(-2*ytmp(initer)*conj(U(initer,initer)));
                additive(initer) = additive(initer) + ytmp(initer);
                for k = initer-1 : -1 : 1
                    additive(k) = additive(k) + ytmp(k);
                    additive = additive - U(:,k)*(2*(U(:,k)'*additive));
                end
                if norm(additive) < eps*norm(x)
                    stag = stag + 1;
                else
                    stag = 0;
                end
                XM = x + additive;
                evalxm = 1;
            elseif evalxm == 1
                addvc = [-(R(1:initer-1,1:initer-1)\R(1:initer-1,initer))*...
                    (w(initer)/R(initer,initer)); w(initer)/R(initer,initer)];
                if norm(addvc) < eps*norm(XM)
                    stag = stag + 1;
                else
                    stag = 0;
                end
                additive = U(:,initer)*(-2*addvc(initer)*conj(U(initer,initer)));
                additive(initer) = additive(initer) + addvc(initer);
                for k = initer-1 : -1 : 1
                    additive(k) = additive(k) + addvc(k);
                    additive = additive - U(:,k)*(2*(U(:,k)'*additive));
                end
                XM = XM + additive;
            end
            r = b - iterapp('mtimes',afun,atype,afcnstr,XM,varargin{:});
            if norm(r) <= tol*n2b
                x = XM;
                flag = 0;
                iter = [outiter, initer];
                

%                 for i = 1:length(XM)/(nx_fine-1)
%                     U_fine_mat(:,i) = [XM((i-1)*(nx_fine-1)+1:(i)*(nx_fine-1))];
%                 end
%         
%                 U_fine_mat_reduce =  U_fine_mat(:,pick);
%                 U0n = U_fine_mat_reduce;
%                 L2NormErrorGMRES(initer+1) = norm(U_fine_mat_reduce' - uF_seq_mat',2);%sqrt(sum(sum(dT*(uF_seq(end-nx_coarse+2:end,:)-U(end-nx_coarse+2:end,:)).^2)));%
%                 LInfNormErrorGMRES(initer+1) = norm(U_fine_mat_reduce' - uF_seq_mat',inf);%max(abs(uF_seq(end-nx_coarse+2:end,:)-U(end-nx_coarse+2:end,:)));%
    
        
                break
            end
            minv_r = r;
            if existM1
                minv_r = iterapp('mldivide',m1fun,m1type,m1fcnstr,r,varargin{:});
                if ~all(isfinite(minv_r))
                    flag = 2;
                    break
                end
            end
            if existM2
                minv_r = iterapp('mldivide',m2fun,m2type,m2fcnstr,minv_r,varargin{:});
                if ~all(isfinite(minv_r))
                    flag = 2;
                    break
                end
            end
            
            normr_act = norm(minv_r);
            resvec((outiter-1)*inner+initer+1) = normr_act/normr0;
            
            if normr_act <= normrmin
                normrmin = normr_act;
                imin = outiter;
                jmin = initer;
                xmin = XM;
                minupdated = 1;
            end
            
%             if normr_act <= tolb
%                 x = xm;
%                 flag = 0;
%                 iter = [outiter, initer];
%                 break
%             else
%                 if stag >= maxstagsteps && moresteps == 0
%                     stag = 0;
%                 end
%                 moresteps = moresteps + 1;
%                 if moresteps >= maxmsteps
%                     if ~warned
%                         warning(message('MATLAB:gmres:tooSmallTolerance'));
%                     end
%                     flag = 3;
%                     iter = [outiter, initer];
%                     break;
%                 end
%             end

%         end
        
        if normr_act <= normrmin
            normrmin = normr_act;
            imin = outiter;
            jmin = initer;
            minupdated = 1;
        end
        
        if stag >= maxstagsteps
            flag = 3;
            break;
        end
        
        
        
        initer
        XM
        
        
        U_fine_mat = [];
        for i = 1:length(XM)/(nx_fine-1)
            U_fine_mat(:,i) = [XM((i-1)*(nx_fine-1)+1:(i)*(nx_fine-1))];
        end
        
        U_fine_mat_reduce =  U_fine_mat(:,pick);
        U0n = U_fine_mat_reduce;
        L2NormErrorGMRES(initer+1) = sqrt(sum(dX*sum(dT*(U_fine_mat_reduce' - uF_seq_mat').^2)));%sqrt(sum(sum(dT*(uF_seq(end-nx_coarse+2:end,:)-U(end-nx_coarse+2:end,:)).^2)));%
        LInfNormErrorGMRES(initer+1) = norm(U_fine_mat_reduce' - uF_seq_mat',inf);%max(abs(uF_seq(end-nx_coarse+2:end,:)-U(end-nx_coarse+2:end,:)));%
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    end         % ends inner loop
    
    evalxm = 0;
    
    if flag ~= 0
        if minupdated
            idx = jmin;
        else
            idx = initer;
        end
        y = R(1:idx,1:idx) \ w(1:idx);
        additive = U(:,idx)*(-2*y(idx)*conj(U(idx,idx)));
        additive(idx) = additive(idx) + y(idx);
        for k = idx-1 : -1 : 1
            additive(k) = additive(k) + y(k);
            additive = additive - U(:,k)*(2*(U(:,k)'*additive));
        end
        x = x + additive;
        xmin = x;
        r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
        minv_r = r;
        if existM1
            minv_r = iterapp('mldivide',m1fun,m1type,m1fcnstr,r,varargin{:});
            if ~all(isfinite(minv_r))
                flag = 2;
                break
            end
        end
        if existM2
            minv_r = iterapp('mldivide',m2fun,m2type,m2fcnstr,minv_r,varargin{:});
            if ~all(isfinite(minv_r))
                flag = 2;
                break
            end
        end
        normr_act = norm(minv_r);
        r = minv_r;
    end
    
    if normr_act <= normrmin
        xmin = x;
        normrmin = normr_act;
        imin = outiter;
        jmin = initer;
    end
    
    if flag == 3
        break;
    end
    if normr_act <= tolb
        flag = 0;
        iter = [outiter, initer];
        break;
    end
    minupdated = 0;
end         % ends outer loop

LInfNormErrorGMRES = LInfNormErrorGMRES';
L2NormErrorGMRES = L2NormErrorGMRES';





% returned solution is that with minimum residual
if flag == 0
    relres = normr_act / n2minv_b;
else
    x = xmin;
    iter = [imin jmin];
    relres = normr_act / n2minv_b;
end

resvec = resvec(1:(outiter-1)*inner+initer+1);
if flag == 2 && initer ~= 0
    resvec(end) = [];
end

% only display a message if the output flag is not used
if nargout < 2
    if restarted
        itermsg(sprintf('gmres(%d)',restart),tol,maxit,[outiter initer],flag,iter,relres);
    else
        itermsg(sprintf('gmres'),tol,maxit,initer,flag,iter(2),relres);
    end
end

function sgn = scalarsign(d)
sgn = sign(d);
if (sgn == 0)
    sgn = 1;
end
