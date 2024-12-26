close all
clear
clc


%% ---------------------------------------------------------------
%  Parareal algorithm for non- linear problem:  
%         y'(t) = y^2(t) , t in (0,T)
%          y(t=0) = -y0
y0=1;
yExactFunction = @(t)(t-1/y0)^(-1);


%% -----------------------------------------------------------------


% initial condition
u0 = -1;

% coarse propagator
N = 10;     %   number of coarse timesteps
T = 1;      %   T = N*Dt  
Dt = T/N;   %   coarse timesteps

xt = linspace(0,T,N+1);

% fine propagator
m = 10;    % number of fine timesteps on each coarse timestep
dt = Dt/m;

K = 20; % number of Parareal iterations

cmap = hsv(K); % color map

% exact solution
 xtrue = linspace(0,T,N*m+1);
 pick = 1:m:N*m+1;
 uF_true = mdl1(u0,T,N*m,'fine');
 
 
%% ------------------------------------------------------------------
% initial coarse integration
U = mdl1(u0,T,N,'coarse');

 figure
 hold on
 plot(xtrue,uF_true,'k',xt,U,'-*m');
 legend('True solution','Parareal solution at k = 0','Location','NW')
 
L2NormError = [ sqrt(sum(Dt*(uF_true(pick)-U).^2))];
LInfNormError = [max(abs(uF_true(pick)-U))];

% iteration loop
for k = 1: K-1
    % k = k + 1
       % parallel loop - predictor
       for n = 1 : N 
           % n = n + 1
           G0 = mdl1(U(:,n),Dt,1,'coarse');
           F0 = mdl1(U(:,n),Dt,m,'fine');
           G(:,n) = G0(:,end);
           F(:,n) = F0(:,end);
       end % end parallel loop
       
       % sequential loop - corrector
       if k < K
           for n = 1: N
               Gn =  mdl1(U(:,n),Dt,1,'coarse');
               U(:, n + 1) = Gn(:,end) - G(:,n) + F(:,n);
           end % sequential loop
       end
       
%        pause
       plot(xt,U,'Color',cmap(k,:),'Marker','o');
       
%        norm(uF_true(pick)-U,'inf')
       L2NormError(k+1) = sqrt(sum(Dt*(uF_true(pick)-U).^2));
       LInfNormError(k+1) = max(abs(uF_true(pick)-U));      
end  % iteration loop

uF_true(pick)
U
L2NormError=L2NormError'
LInfNormError=LInfNormError'

% figure
% semilogy(0:K-1,LInfNormError,'bx-',0:K-1,L2NormError,'rx-')
% legend('L^{oo}NormError','L^2NormError')
% xlabel('k')

figure
semilogy(0:K-1,LInfNormError,'b--^')
legend('L^{\infty}NormError')
xlabel('k')
title(['T=',num2str(T)])
