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
N = 20;     %   number of coarse timesteps
T = 1;      %   T = N*Dt  
Dt = T/N;   %   coarse timesteps

xt = linspace(0,T,N+1);

% fine propagator
m = 2;    % number of fine timesteps on each coarse timestep
dt = Dt/m;

K = 10 % number of Parareal iterations

cmap = hsv(K); % color map

% exact solution
 xtrue = linspace(0,T,N*m+1);
 pick = 1:m:N*m+1;
 uF_true = mdl(u0,T,N*m,'fine');
 

 
%% ------------------------------------------------------------------
% initial coarse integration
U = mdl(u0,T,N,'coarse');

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
FC_boundL2 = L2NormError;

phi = 1/(1-a*dt);
phiDT = 1/(1-a*Dt);
% iteration loop
for k = 1: K-1
    % k = k + 1
       % parallel loop - predictor
       for n = 1 : N 
           % n = n + 1
           G0 = mdl(U(:,n),Dt,1,'coarse');
           F0 = mdl(U(:,n),Dt,m,'fine');
           G(:,n) = G0(:,end);
           F(:,n) = F0(:,end);
       end % end parallel loop
       
       % sequential loop - corrector
       if k < K
           for n = 1: N
               Gn =  mdl(U(:,n),Dt,1,'coarse');
               U(:, n + 1) = Gn(:,end) - G(:,n) + F(:,n);
           end % sequential loop
       end
        
%        pause
       plot(xt,U,'Color',cmap(k,:),'Marker','o');

%        
%        norm(uF_true(pick)-U,'inf')
       proNj(k+1) = proNj(k)*(N-k);
       L2NormError(k+1) = sqrt(sum(Dt*(uF_true(pick)-U).^2));
       LInfNormError(k+1) = max(abs(uF_true(pick)-U));
       SupperlinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt)))^(k)/factorial(k)*proNj(k+1)*max(abs(uF_true(pick)-U0n));
       LinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt))/(1- abs(Rz(a*Dt))))^(k)*max(abs(uF_true(pick)-U0n));
     % FC bound L2 
      FC_boundL2(k+1) = (abs(phi^m - phiDT)*(1 - abs(phiDT)^((N-1)/m) )/ (1 - abs(phiDT)))^(k)*sqrt(sum(Dt*(uF_true(pick)-U0n).^2));

end  % iteration loop

uF_true_coarse=uF_true(pick)
U

L2NormError=L2NormError'
LInfNormError=LInfNormError'
SupperlinearErrorBound=SupperlinearErrorBound'
LinearErrorBound=LinearErrorBound'

% figure
% semilogy(0:K-1,LInfNormError,'b--^',0:K-1,L2NormError,'g--^',0:K-1,SupperlinearErrorBound,'r*-',0:K-1,LinearErrorBound,'mx-')
% legend('L^{oo}NormError','L^2NormError','Superlinear bound','Linear bound')
% title(['T=',num2str(T)])
% xlabel('k')
% ylim([1e-18 1e0])
% %FC bound L2

figure
semilogy(0:K-1,L2NormError,'b--^',0:K-1,SupperlinearErrorBound,'k.-',0:K-1,FC_boundL2,'ro-')
legend('L^2- Parareal error norm','Superlinear bound','FC bound L^2- norm')
xlabel('Iteration ')
ylabel('Error')
ylim([1e-18 1])
xlim([0 7])
title(['T=',num2str(T)])
