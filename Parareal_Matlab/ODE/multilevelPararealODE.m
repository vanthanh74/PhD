close all
clear
clc


%% ---------------------------------------------------------------
%  Parareal algorithm for linear problem:  
%         y'(t) = a*y(t), t in (0,T)
%          y(t=0) = y0

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
N = 10;     %   number of coarse timesteps
T = 1;      %   T = N*Dt  
Dt = T/N;   %   coarse timesteps

xt = linspace(0,T,N+1);

% fine propagator
m = 20;    % number of fine time steps on each coarse time interval
dt = Dt/m;

K = 7 % number of Parareal iterations

cmap = hsv(K); % color map
%%
% Reduced system
F_tilde = (1-a*dt)^-m;  % fine propagator
G_tilde = (1-a*Dt)^-1;  % coarse propagator

Uk_0 = zeros(N,1);
rhs0 = [u0; zeros(m*N,1)];
rhs = [u0; (F_tilde-G_tilde)*Uk_0];

A = eye(m*N+1,m*N+1) + (-(1-a*dt)^-1)*diag(ones(m*N,1),-1);  % original problem matrix

A_tilde = eye(N+1,N+1) + (-F_tilde)*diag(ones(N,1),-1);      % reduced matrix

M_tilde = eye(N+1,N+1) + (-G_tilde)*diag(ones(N,1),-1);      % approximate reduced matrix

% exact solution

 xtrue = linspace(0,T,N*m+1);
 pick = 1:m:N*m+1;
 
 uF_true = A\rhs0;
 
  
%% ------------------------------------------------------------------
% initial coarse integration

U = M_tilde\rhs;

 figure
 hold on
 plot(xtrue,uF_true,'k',xt,U,'-*m');
 legend('True solution','Parareal solution at k = 0')
 
U0n = U;
proNj = 1;

L2NormError = [ sqrt(sum(Dt*(uF_true(pick)-U).^2))];
LInfNormError = [max(abs(uF_true(pick)-U))];
SupperlinearErrorBound  = LInfNormError;
LinearErrorBound  = LInfNormError;

%%



% iteration loop
for k = 1: K-1
        % k = k + 1
        rhs(2:end) =  [ (F_tilde-G_tilde)*U(1:end-1)];
        U = M_tilde\rhs;

       plot(xt,U,'Color',cmap(k,:),'Marker','o');
       
%        
%        norm(uF_true(pick)-U,'inf')
       proNj(k+1) = proNj(k)*(N-k);
       L2NormError(k+1) = sqrt(sum(Dt*(uF_true(pick)-U).^2));
       LInfNormError(k+1) = max(abs(uF_true(pick)-U));
       SupperlinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt)))^(k)/factorial(k)*proNj(k+1)*max(abs(uF_true(pick)-U0n));
       LinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt))/(1- abs(Rz(a*Dt))))^(k)*max(abs(uF_true(pick)-U0n));
end  % iteration loop

uF_true_coarse=uF_true(pick)
U

L2NormError=L2NormError'
LInfNormError=LInfNormError'
SupperlinearErrorBound=SupperlinearErrorBound'
LinearErrorBound=LinearErrorBound'

figure
semilogy(0:K-1,LInfNormError,'b--',0:K-1,L2NormError,'g--',0:K-1,SupperlinearErrorBound,'r*-',0:K-1,LinearErrorBound,'mx-')
legend('L^{oo}NormError','L^2NormError','Superlinear bound','Linear bound')
xlabel('k')





















