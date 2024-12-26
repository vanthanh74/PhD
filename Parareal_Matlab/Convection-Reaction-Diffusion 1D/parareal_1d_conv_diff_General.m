function parareal_1d_conv_diff_General(L,T,dT,dX,dt,dx,nt_interval,nt_coarse,nx_coarse,nt_fine,nx_fine,M,m,K,solver)


% Illustration of the Parareal method (using the backward Euler
% discretization for coarse and fine propagator).

% IVP: u_t = a*u_xx - b*u_x + c*u + f               in (0,L) x (0,T),
%      u(x,0) = u0(x)                    in (0,L),
%      u(0,t) = u(L,t) = 0                 t in (0,T)

% INPUT PARAMETERS:
%   * a>0: constant problem parameter
%   * nt_coarse: number of coarse intervals in time
%   * nx_coarse: number of coarse intervals in space
%   * nt_fine: number of fine intervals per coarse interval in time
%   * nx_fine: number of fine intervals per coarse interval in space
%   * K: number of Parareal iterations (2 <= K <= n_coarse+1)
%        (--> K sequential coarse sweeps and K-1 parallel fine sweeps)
% OUTPUT:
%   * LInfinityErrorNorm
%
% Organization of memory:
% * Coarse solution: vector of length n_coarse+1
% * Fine reference solution: vector of length n_coarse*n_fine+1
% * Parallel iterates: (ncoarse x n_fine+1)-array

global xC xF a b c

% For the plots


% coarse points
xC = linspace(0,L,nx_coarse+1);
tC = linspace(0,T,nt_interval+1);

% dX = L/nx_coarse; % coarse spatial discretization steps
% dT = T/nt_coarse; % coarse temporal discretization steps

% fine points
xF = linspace(0,L,nx_fine+1);
tF = linspace(0,T,nt_fine+1);


cmap = hsv(K);

% %--------------------------------------
% exact solution 1
u_ex=@(x,t) x*(L-x)^2*exp(-2*t);

% Initial condition u(x,0) = u0(x)
ux0_function = @(x) (x.*(L-x).^2);%x.^4.*(1-x);%;

%     ux0_function = @(x) zeros(1,length(x));

%syms a b c x t L
% u = x*(L-x)^2*exp(-2*t)
% f = diff(u,t) - a*diff(diff(u,x),x) + b*diff(u,x) - c*u
% f = @(x,t) 2*exp(-2*t)*(2*L - 2*x) - 2*x*exp(-2*t) - 2*x*exp(-2*t)*(L - x)^2;%x^4*(1-x) + t^2;%

f = @(x,t) b*(exp(-2*t)*(L - x)^2 - x*exp(-2*t)*(2*L - 2*x)) - a*(2*x*exp(-2*t) - 2*exp(-2*t)*(2*L - 2*x)) - 2*x*exp(-2*t)*(L - x)^2 - c*x*exp(-2*t)*(L - x)^2;
% initial solution u0
ux0 = ux0_function(xF);

% % %--------------------------------------
% % exact solution 2
% u_ex=@(x,t) sin(2*pi*x)*exp(-2*t);
% %syms a b c x t L
% % u = sin(2*pi*x)*exp(-2*t);%x*(L-x)^2*exp(-2*t)
% % f = diff(u,t) - a*diff(diff(u,x),x) + b*diff(u,x) - c*u
% ux0_function = @(x) sin(2*pi*x);
% f = @(x,t) 4*a*pi^2*exp(-2*t)*sin(2*pi*x) - c*exp(-2*t)*sin(2*pi*x) - 2*exp(-2*t)*sin(2*pi*x) + 2*b*pi*exp(-2*t)*cos(2*pi*x);
% ux0 = ux0_function(xF);













%                           ux0
%       ----*----*----*----*----*----*----*---->x
%      |
%      |
%      |
%      |
%  u0t |
%      |
%      |
%      |
%      |
%      t





% if strcmp(solver,'B')==1

%% ------------------------------------------------------------------
% initial coarse integration
ux0_restriction = restriction(ux0,nx_coarse,nx_fine);
U = [ux0_restriction];
U1=U;
U_interpolation = [ux0];
U_interpolation1=U_interpolation;
for i =1:nt_interval
    %      i = i+1
    U_temp = BackwardEuler_convdiff(xC,dX,dT,f,U(i,:),tC(i),nx_coarse,M);
    U(i+1,:) = U_temp(end,:);
    U_interpolation(i+1,:) = spline(xC,U(i+1,:),xF);
end

% figure
% surf(U)
% figure
% surf(U_interpolation)

%     figure
%     hold on
%     plot(xF,U_interpolation(1,:))
%
% Sequential fine solution
U_FineSequential=[ux0];
for n = 1:nt_interval
    if strcmp(solver,'BackwardEuler')==1
        U_Fine = BackwardEuler_convdiff(xF,dx,dt,f,U_FineSequential(n,:),tC(n),nx_fine,m) ;
    elseif strcmp(solver,'Runge-Kutta4')==1
        U_Fine = RungeKuttaFine_convdiff(dx,dt,f,U_FineSequential(n,:),tC(n),nx_fine,m) ;
    end
    U_FineSequential(n+1,:) = U_Fine(end,:);
    %             plot(xF,U_FineSequential(n,:))
    %             pause(0.5)
end

% Exact solution at time T = 1
U_exact=zeros(1,nx_fine+1);
for i=1:nx_fine+1
    U_exact(i)=u_ex(xF(i),T);
end

% norm(U_interpolation(end,:) - U_exact,inf);


% LInfNormError = zeros(K,1);
% LInfNormError(1) = norm(uF_true(t_pick,x_pick)-U_interpolation,inf);

normInf2executiveParearealIter = [];
normInfDefect = [norm(U_interpolation - U_FineSequential,inf)];%[max(max(abs(U_interpolation - U_FineSequential)))];%

L2NormError = norm(U_FineSequential-U_interpolation,2);%[ sqrt(sum(sum((U_FineSequential-U_interpolation).^2)))];%;%norm(U_FineSequential(end,:)-U_interpolation(end,:),2);%
LInfNormError  = norm(U_interpolation-U_FineSequential,inf);%[max(max(abs(U_FineSequential-U_interpolation)))];%

U0 = U_interpolation;

cmap = hsv(K); % color map

% SupperlinearErrorBound  = normInfDefect;
% LinearErrorBound  = normInfDefect;
% proNj = 1;

figure
hold on
plot(xF,U_interpolation(end,:),'-*m',xF,U_exact,'-*k');
legend('U0','U_{exact}')
title(['Solution at time T = ',num2str(T)])
xlabel('x')
ylabel('u')



alpha = 1;
% iteration loop
for k = 1: K-1
    % k = k + 1
    % parallel loop - predictor
    %         parpool(nt_coarse)
    for n = 1 : nt_interval
        % n = n + 1
        if strcmp(solver,'BackwardEuler')==1
            F0 =   BackwardEuler_convdiff(xF,dx,dt,f,U_interpolation(n,:),tC(n),nx_fine,m) ;
        elseif strcmp(solver,'Runge-Kutta4')==1
            F0 =   RungeKuttaFine_convdiff(dx,dt,f,U_interpolation(n,:),tC(n),nx_fine,m) ;
        end
        ux0F_restriction = restriction(U_interpolation(n,:),nx_coarse,nx_fine);
        G0 =  BackwardEuler_convdiff(xC,dX,dT,f,ux0F_restriction,tC(n),nx_coarse,M) ;
        G(n,:) = G0(end,:);
        F(n,:) = F0(end,:);
        
    end % end parallel loop
    
    
    % sequential loop - corrector
    if k < K
        for n = 1: nt_interval
            %                F0 =  first_order_upwind_adv(a,U(n,:),u0t(n+1),dt,dx,m,'fine');
            ux0F_restriction = restriction(U_interpolation(n,:),nx_coarse,nx_fine);
            Gn = BackwardEuler_convdiff(xC,dX,dT,f,ux0F_restriction,tC(n),nx_coarse,M) ;
            U_interpolation( n + 1,:) = alpha*(  spline(xC,Gn(end,:),xF)- spline(xC,G(n,:),xF)  )  + F(n,:)     ;%...
            %                                             + 0.05*( F(n+1,:) - spline(xC,G(n+1,:),xF));
            %                U( n + 1,:) = restriction(U_temp,nx_coarse);
            norm((spline(xC,Gn(end,:),xF)  - spline(xC,G(n,:),xF)));
        end % sequential loop
    end
    
    L2NormError(k+1)    = norm(U_FineSequential-U_interpolation,2);%[ sqrt(sum(sum((U_FineSequential-U_interpolation).^2)))];%
    LInfNormError(k+1)  = norm(U_FineSequential-U_interpolation,inf);%[max(max(abs(U_FineSequential-U_interpolation)))];%
    
    normInf_U0_Uinterpolate = norm(U0-U_interpolation,inf);%max(max(abs(U0-U_interpolation))) ;%norm(U0(end,:)-U_interpolation(end,:),inf) ;%
    normInf2executiveParearealIter =[normInf2executiveParearealIter;normInf_U0_Uinterpolate] ;
    normInfDefect =[normInfDefect ; norm(U_interpolation-U_FineSequential,inf)];%[normInfDefect ; max(max(abs(U_interpolation-U_FineSequential)))];%[normInfDefect ; norm(U_interpolation(end,:)-U_FineSequential(end,:),inf)];%
    
    %     SupperlinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt)))^(k)/factorial(k)*proNj(k+1)*max(abs(uF_true(pick)-U0n));
    %        LinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt))/(1- abs(Rz(a*Dt))))^(k)*max(abs(uF_true(pick)-U0n));
    
    %       figure
    %       surf(xF,tC,U_interpolation)
    
    U0 = U_interpolation;
    %         pause(0.2)
    plot(xF,U_interpolation(end,:),'Color',cmap(k,:),'Marker','o');
    %      proNj(k+1) = proNj(k)*(N-k);.
    
    
end  % iteration loop
%     normInf2executiveParearealIter
%     normInfDefect
%     LInfNormError = 1;

figure
surf(xF,tC,U_interpolation)
xlabel('x')
ylabel('t')
zlabel('u')
title('Solution')


%     figure
%     surf(xF,tC,U_FineSequential)
%     xlabel('x')
%     ylabel('t')
%     zlabel('u')
%     title('Fine Sequential Solution')
%      norm(  U_FineSequential -  U_interpolation)

figure
semilogy(1:K-1, normInf2executiveParearealIter,'b--^',0:K-1,normInfDefect,'r--^')
legend('normL^{\infty}2consecutiveParearealIter','normL^{\infty}Defect')
xlabel('k')
title(solver)


figure
semilogy(0:K-1,LInfNormError,'r-^',0:K-1,L2NormError,'b-o')
legend('L^{\infty}NormError','L^2NormError')
xlabel('k')
title('Error norm between the Parareal solution and the Sequential Fine solution at the end of each time-slice')
% % ylim([1e-16 1])

% figure
% semilogy(0:K-1,LInfNormError,'b-o')
% legend('L^{\infty}NormError')
% xlabel('k')
% title('Error norm between the Parareal solution and the Sequential Fine solution at the end of each time-slice')
% ylim([1e-16 1])
