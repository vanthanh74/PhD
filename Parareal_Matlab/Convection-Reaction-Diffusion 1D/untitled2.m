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

T = 1; %  Intervall (0,T)
L = 1; % omega=(0,L)

a = 3;
b = 0.0005;
c = 1;

nt_interval = 4; % number of time intervals

nx_coarse = 5; % number of coarse spatial intervals
nt_coarse = 4;% 300 %20    % number of coarse time intervals


nx_fine = 	5; % number of fine spatial intervals
nt_fine = 8;% 3000 %200 . % number of fine time intervals


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

K = 7    % number of Parareal iterations


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

% exact solution
u_ex=@(x,t) x*(L-x)^2*exp(-2*t);

ux0_function = @(x) (x.*(L-x).^2);

f = @(x,t) b.*(exp(-2.*t).*(L - x).^2 - x.*exp(-2.*t).*(2.*L - 2.*x)) - a.*(2.*x.*exp(-2.*t) - 2.*exp(-2.*t).*(2.*L - 2.*x)) - 2.*x.*exp(-2.*t).*(L - x).^2 - c.*x.*exp(-2.*t).*(L - x).^2;



%% Define find and coarse propagators
% Backward Euler fine solver
% A_tilde
rF = a*dt/(dx^2);
kF = b*dt/(2*dx);


A_F = (kF-rF)*diag(ones(nx_fine-2,1),1) + (1+2*rF-c*dt)*eye(nx_fine-1) - (kF+rF)*diag(ones(nx_fine-2,1),-1);
Id_F = eye(size(A_F));
A_F_Global =  zeros((nt_interval+1)*(nx_fine-1));
for i =1:nt_interval+1
    %  i = i + 1
    A_F_Global((i-1)*(nx_fine-1)+1:i*(nx_fine-1),(i-1)*(nx_fine-1)+1:i*(nx_fine-1)) = Id_F;
    if i>1
        A_F_Global((i-1)*(nx_fine-1)+1:i*(nx_fine-1),(i-2)*(nx_fine-1)+1:(i-1)*(nx_fine-1)) = -A_F^-m;
    end
end
A_tilde = A_F_Global;

% Backward Euler coarse solver
% M_tilde
rC = a*dT/(dX^2);
kC = b*dT/(2*dX);

A_C = (kC-rC)*diag(ones(nx_coarse-2,1),1) + (1+2*rC-c*dT)*eye(nx_coarse-1) - (kC+rC)*diag(ones(nx_coarse-2,1),-1);
Id_C = eye(size(A_C));
A_C_Global =  zeros((nt_interval+1)*(nx_coarse-1));
for i =1:nt_coarse+1
    %  i = i + 1
    A_C_Global((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i-1)*(nx_coarse-1)+1:i*(nx_coarse-1)) = Id_C;
    if i>1
        A_C_Global((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i-2)*(nx_coarse-1)+1:(i-1)*(nx_coarse-1)) = -A_C^-1;
    end
end
M_tilde = A_C_Global;

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
for i=1:nx_fine+1
    U_exact(i)=u_ex(xF(i),T);
end

% initial solution u0
ux0 = ux0_function(xF(2:end-1))';

% Fine sequential propagator
Uk_0 = zeros((nt_interval+1)*(nx_fine-1),1);
Uk_0(1:nx_fine-1) = ux0;
rhs0 = [ux0; zeros(((nt_interval+1)*(nx_fine-1)-nx_fine+1),1)];
%  rhs = [ux0; (A_tilde-M_tilde)*Uk_0];
% source force F
F = zeros((nt_interval+1)*(nx_fine-1),1);
F(1:nx_fine-1) = 0;
 i_temp = 1;
t0 = tF(1);
% for j = 1: nt_fine
%     t0 = t0 + dt;
%     i_temp = 1;
    for i = 2:nt_interval+1
        F(i_temp*(nx_fine-1) + 1 :i_temp*(nx_fine-1) + nx_fine-1,:) = A_F^-m*dt*( f(xF(2:nx_fine),tF(i))' +  A_F*f(xF(2:nx_fine),tC(i))' );
        i_temp = i_temp + 1;
    end
    uF_true = A_tilde\(Uk_0 + F);
%     Uk_0 = uF_true;
% end
figure
hold on
plot(xF,[0;uF_true(end-nx_fine+2:end);0])

plot(xF,U_exact)
legend('U_{FineSequential}','U_{ex}')
title('Sequentially Fine Solution and the Exact Solution')



%% ------------------------------------------------------------------
% initial coarse integration

% initial solution u0
ux0C = ux0_function(xC(2:end-1))';

% Coarse sequential propagator

Uk_0C = zeros((nt_coarse+1)*(nx_coarse-1),1);
Uk_0C(1:nx_coarse-1) = ux0C;
rhs0 = [ux0; zeros(((nt_coarse+1)*(nx_coarse-1)-nx_coarse+1),1)];
%  rhs = [ux0; (A_tilde-M_tilde)*Uk_0];
% source force F
F_C = zeros((nt_coarse+1)*(nx_coarse-1),1);
F_C(1:nx_coarse-1) = 0;
i_temp = 1;
for i = 2:nt_coarse+1
    F_C(i_temp*(nx_coarse-1) + 1 :i_temp*(nx_coarse-1) +  nx_coarse-1,:) = A_C^-1*dT*f(xC(2:nx_coarse),tC(i))';
    i_temp = i_temp + 1;
end

figure
hold on
uC = M_tilde\(Uk_0C + F_C);
plot(xF,[0;uC(end-nx_coarse+2:end);0])

plot(xF,U_exact)
legend('U_{0}','U_{ex}')
title('Parareal Solution and the Exact Solution')

U0n = uC;
% proNj = 1;

L2NormError = [ sqrt(sum(dT*(uF_true-U).^2))];
LInfNormError = [max(abs(uF_true-U))];
% SupperlinearErrorBound  = LInfNormError;
% LinearErrorBound  = LInfNormError;

%%


% iteration loop
for k = 1: K-1
    % k = k + 1
    %         rhs(2:end) =  [ (F_tilde-G_tilde)*U(1:end-1)];
    %         U = M_tilde\rhs;
    U = M_tilde\((M_tilde - A_tilde)*uC + rhs);
    
    plot(xt,U,'Color',cmap(k,:),'Marker','o');
    
    %
    %        norm(uF_true(pick)-U,'inf')
%     proNj(k+1) = proNj(k)*(N-k);
    L2NormError(k+1) = sqrt(sum(Dt*(uF_true-U).^2));
    LInfNormError(k+1) = max(abs(uF_true-U));
%     SupperlinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt)))^(k)/factorial(k)*proNj(k+1)*max(abs(uF_true-U0n));
%     LinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt))/(1- abs(Rz(a*Dt))))^(k)*max(abs(uF_true-U0n));
end  % iteration loop

uF_true_coarse=uF_true(pick)
U

L2NormError=L2NormError'
LInfNormError=LInfNormError'
SupperlinearErrorBound=SupperlinearErrorBound'
LinearErrorBound=LinearErrorBound'

figure
semilogy(0:K-1,LInfNormError,'b--^',0:K-1,L2NormError,'g--^',0:K-1,SupperlinearErrorBound,'r*-',0:K-1,LinearErrorBound,'mx-')
legend('L^{oo}NormError','L^2NormError','Superlinear bound','Linear bound')
xlabel('k')
ylim([1e-18 1])




















