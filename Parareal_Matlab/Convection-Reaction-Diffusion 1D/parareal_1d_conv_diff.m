function parareal_1d_conv_diff(L,T,dT,dX,dt,dx,nt_coarse,nx_coarse,nt_fine,nx_fine,m,K,solver)

% Illustration of the Parareal method (using the backward Euler
% discretization for coarse and fine propagator).

% IVP: u_t = u_xx - b*u_x + c*u + f                   in (0,L) x (0,T),
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
tC = linspace(0,T,nt_coarse+1);

% dX = L/nx_coarse; % coarse spatial discretization steps
% dT = T/nt_coarse; % coarse temporal discretization steps

% fine points
xF = linspace(0,L,nx_fine+1);
tF = linspace(0,T,nt_fine+1);


cmap = hsv(K);

% exact solution
u_ex=@(x,t) x*(L-x)^2*exp(-2*t);
% Initial condition u(x,0) = u0(x)
ux0_function = @(x) (x.*(L-x).^2);%x.^4.*(1-x);%;
%     ux0_function = @(x) zeros(1,length(x));

%syms a b c x t L
% u = x*(L-x)^2*exp(-2*t)
% f = diff(u,t) - a*diff(diff(u,x),x) + b*diff(u,x) - c*u
% f = @(x,t) 2*exp(-2*t)*(2*L - 2*x) - 2*x*exp(-2*t) - 2*x*exp(-2*t)*(L - x)^2;%x^4*(1-x) + t^2;%
f = @(x,t) b*(exp(-2*t)*(L - x)^2 - x*exp(-2*t)*(2*L - 2*x)) - a*(2*x*exp(-2*t) - 2*exp(-2*t)*(2*L - 2*x)) - 2*x*exp(-2*t)*(L - x)^2 - c*x*exp(-2*t)*(L - x)^2;
% f = @(x,t) x^4*(1-x) +t^2;
% f = @(x,t) 2*a*exp(-2*t) - b*(x*exp(-2*t) - exp(-2*t)*(L - x)) - 2*x*exp(-2*t)*(L - x);
% initial solution u0
ux0 = ux0_function(xF);
%     u0t = u0t_function(tF)';


%     u0ttrue = u0t_function(ttrue)';


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
    % initial fine integration
    ux0_restriction = restriction(ux0,nx_coarse,nx_fine);
    U = [ux0_restriction];
    U_interpolation = [ux0];
    for i =1:nt_coarse
        %      i = i+1
        U(i+1,:) = BackwardEulerCoarse_convdiff(dX,dT,f,U(i,:),tC(i),nx_coarse);
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
    for n = 1:nt_coarse
          if strcmp(solver,'BackwardEuler')==1
               U_Fine = BackwardEulerFine_convdiff(dx,dt,f,U_FineSequential(n,:),tC(n),nx_fine,m) ;
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
    % L2NormError = [ sqrt(sum(dT*(uF_true(pick)-U).^2))];
    % LInfNormError = [max(abs(uF_true(pick)-U))];
    % LInfNormError = zeros(K,1);
    % LInfNormError(1) = norm(uF_true(t_pick,x_pick)-U_interpolation,inf);
    normInf2executiveParearealIter = [];
    normInfDefect = [norm(U_interpolation(end,:) - U_FineSequential(end,:),inf)];
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
        for n = 1 : nt_coarse
            % n = n + 1
            if strcmp(solver,'BackwardEuler')==1
                 F0 =   BackwardEulerFine_convdiff(dx,dt,f,U_interpolation(n,:),tC(n),nx_fine,m) ;
            elseif strcmp(solver,'Runge-Kutta4')==1
                 F0 =   RungeKuttaFine_convdiff(dx,dt,f,U_interpolation(n,:),tC(n),nx_fine,m) ;
            end
            ux0F_restriction = restriction(U_interpolation(n,:),nx_coarse,nx_fine);
            G0 =  BackwardEulerCoarse_convdiff(dX,dT,f,ux0F_restriction,tC(n),nx_coarse)' ;
            G(n,:) = G0(end,:);
            F(n,:) = F0(end,:);

        end % end parallel loop
         
           
        % sequential loop - corrector
        if k < K
            for n = 1: nt_coarse
                %                F0 =  first_order_upwind_adv(a,U(n,:),u0t(n+1),dt,dx,m,'fine');
                ux0F_restriction = restriction(U_interpolation(n,:),nx_coarse,nx_fine);
                Gn = BackwardEulerCoarse_convdiff(dX,dT,f,ux0F_restriction,tC(n),nx_coarse)' ;
                U_interpolation( n + 1,:) = alpha*(  spline(xC,Gn(end,:),xF)- spline(xC,G(n,:),xF)  )  + F(n,:)     ;%...
%                                             + 0.05*( F(n+1,:) - spline(xC,G(n+1,:),xF));
                %                U( n + 1,:) = restriction(U_temp,nx_coarse);
                norm((spline(xC,Gn(end,:),xF)  - spline(xC,G(n,:),xF)));
            end % sequential loop
        end
        normInf_U0_Uinterpolate = norm(U0(end,:)-U_interpolation(end,:),inf) ;
        normInf2executiveParearealIter =[normInf2executiveParearealIter;normInf_U0_Uinterpolate] ;
        normInfDefect =[normInfDefect ; norm(U_interpolation(end,:)-U_FineSequential(end,:),inf)];
        %     SupperlinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt)))^(k)/factorial(k)*proNj(k+1)*max(abs(uF_true(pick)-U0n));
        %        LinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt))/(1- abs(Rz(a*Dt))))^(k)*max(abs(uF_true(pick)-U0n));
        
        %       figure
        %       surf(xF,tC,U_interpolation)
        
        U0 = U_interpolation;
        plot(xF,U_interpolation(end,:),'Color',cmap(k,:),'Marker','o');
        %      proNj(k+1) = proNj(k)*(N-k);
        
        
    end  % iteration loop
    normInf2executiveParearealIter
    normInfDefect
    LInfNormError = 1;
    
    figure
    surf(xF,tC,U_interpolation)
    xlabel('x')
    ylabel('t')
    zlabel('u')
    

    figure
    semilogy(1:K-1, normInf2executiveParearealIter,'-*b',0:K-1,normInfDefect,'-*r')
    legend('normInf2executiveParearealIter','normInfDefect')
    xlabel('k')
    title(solver)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
% elseif strcmp(solver,'F')==1
%     
%     
%     
%     
%     %% ------------------------------------------------------------------
%     % initial fine integration
%     ux0_restriction = restriction(ux0,nx_coarse,nx_fine);
%     U = [ux0_restriction];
%     U_interpolation = [ux0];
%     for i =1:nt_coarse
%         %      i = i+1
%         U(i+1,:) = ForwardEulerCoarse_heat(dX,dT,f,U(i,:),tC(i),nx_coarse);
%         U_interpolation(i+1,:) = spline(xC,U(i+1,:),xF);
%     end
%     
%     % figure
%     % surf(U)
%     % figure
%     % surf(U_interpolation)
%     
%     % Sequential fine solution
%     U_FineSequential=[ux0];
%     for n = 1:nt_coarse
%         U_Fine = ForwardEulerFine_heat(dx,dt,f,U_FineSequential(n,:),tC(n),nx_fine,m) ;
%         U_FineSequential(n+1,:) = U_Fine(end,:);
%     end
%     
%     % Exact solution at time T = 1
%     U_exact=zeros(1,nx_fine+1);
%     for i=1:nx_fine+1
%         U_exact(i)=u_ex(xF(i),T);
%     end
%     
%     % norm(U_interpolation(end,:) - U_exact,inf);
%     % L2NormError = [ sqrt(sum(dT*(uF_true(pick)-U).^2))];
%     % LInfNormError = [max(abs(uF_true(pick)-U))];
%     % LInfNormError = zeros(K,1);
%     % LInfNormError(1) = norm(uF_true(t_pick,x_pick)-U_interpolation,inf);
%     normInf2executedParearealIter = [];
%     normInfDefect = [norm(U_interpolation(end,:) - U_FineSequential(end,:),inf)];
%     U0 = U_interpolation;
%     
%     cmap = hsv(K); % color map
%     
%     % SupperlinearErrorBound  = normInfDefect;
%     % LinearErrorBound  = normInfDefect;
%     % proNj = 1;
%     
%     figure
%     hold on
%     plot(xF,U_interpolation(end,:),'-*m',xF,U_exact,'-*k');
%     legend('U0','U_{exact}')
%     title(['Solution at time T = ',num2str(T)])
%     xlabel('x')
%     ylabel('u')
%     % iteration loop
%     for k = 1: K-1
%         % k = k + 1
%         % parallel loop - predictor
%         for n = 1 : nt_coarse
%             % n = n + 1
%             F0 =   ForwardEulerFine_heat(dx,dt,f,U_interpolation(n,:),tC(n),nx_fine,m) ;
%             ux0F_restriction = restriction(U_interpolation(n,:),nx_coarse,nx_fine);
%             G0 =  ForwardEulerCoarse_heat(dX,dT,f,ux0F_restriction,tC(n),nx_coarse)' ;
%             G(n,:) = G0(end,:);
%             F(n,:) = F0(end,:);
%         end % end parallel loop
%         
%         
%         
%         % sequential loop - corrector
%         if k < K
%             for n = 1: nt_coarse
%                 %                F0 =  first_order_upwind_adv(a,U(n,:),u0t(n+1),dt,dx,m,'fine');
%                 ux0F_restriction = restriction(U_interpolation(n,:),nx_coarse,nx_fine);
%                 Gn = ForwardEulerCoarse_heat(dX,dT,f,ux0F_restriction,tC(n),nx_coarse)' ;
%                 U_interpolation( n + 1,:) = spline(xC,Gn(end,:),xF)  - spline(xC,G(n,:),xF) + F(n,:);
%                 %                U( n + 1,:) = restriction(U_temp,nx_coarse);
%             end % sequential loop
%         end
%         normInf_U0_Uinterpolate = norm(U0(end,:)-U_interpolation(end,:),inf) ;
%         normInf2executedParearealIter =[normInf2executedParearealIter;normInf_U0_Uinterpolate] ;
%         normInfDefect =[normInfDefect ; norm(U_interpolation(end,:)-U_FineSequential(end,:),inf)];
%         %     SupperlinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt)))^(k)/factorial(k)*proNj(k+1)*max(abs(uF_true(pick)-U0n));
%         %        LinearErrorBound(k+1)  = (abs(exp(a*Dt) - Rz(a*Dt))/(1- abs(Rz(a*Dt))))^(k)*max(abs(uF_true(pick)-U0n));
%         
%         %       figure
%         %       surf(xF,tC,U_interpolation)
%         
%         U0 = U_interpolation;
%         plot(xF,U_interpolation(end,:),'Color',cmap(k,:),'Marker','o');
%         %      proNj(k+1) = proNj(k)*(N-k);
%         
%         
%     end  % iteration loop
%     normInf2executedParearealIter
%     normInfDefect
%    
%     
%     figure
%     semilogy(1:K-1, normInf2executedParearealIter,'-*b',0:K-1,normInfDefect,'-*r')
%     legend('normInf2executedParearealIter','normInfDefect')
%     xlabel('k')
%     title('ForwardEuler')
% end



