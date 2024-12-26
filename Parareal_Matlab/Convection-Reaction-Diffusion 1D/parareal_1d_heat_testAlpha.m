function   parareal_1d_heat_testAlpha(L,T,dT,dX,dt,dx,nt_coarse,nx_coarse,nt_fine,nx_fine,m,K)

% Illustration of the Parareal method (using the backward Euler
% discretization for coarse and fine propagator).

% IVP: u_t = u_xx + f               in (0,L) x (0,T),
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

global xC xF

% For the plots


% coarse points
xC = linspace(0,L,nx_coarse+1);
tC = linspace(0,T,nt_coarse+1);

% dX = L/nx_coarse; % coarse spatial discretization steps
% dT = T/nt_coarse; % coarse temporal discretization steps

% fine points
xF = linspace(0,L,nx_fine+1);
tF = linspace(0,T,nt_fine+1);

xtrue = linspace(0,L,nx_coarse*nx_fine+1);
ttrue = linspace(0,T,nt_coarse*nt_fine+1);

% dx = L/nx_fine;   % fine spatial discretization steps
% dt = dT/m;   % fine temporal discretization steps



cmap = hsv(K);

% exact solution
u_ex=@(x,t) x*(L-x)^2*exp(-2*t);
% Initial condition u(x,0) = u0(x)
ux0_function = @(x) (x.*(L-x).^2);
%     u0t_function = @(t) 0;
f = @(x,t) 2*exp(-2*t)*(2*L - 2*x) - 2*x*exp(-2*t) - 2*x*exp(-2*t)*(L - x)^2;
% initial solution u0
ux0 = ux0_function(xF);
%     u0t = u0t_function(tF)';

ux0true = ux0_function(xtrue);
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






    
    %% ------------------------------------------------------------------
    % initial fine integration
    ux0_restriction = restriction(ux0,nx_coarse,nx_fine);
    U = [ux0_restriction];
    U_interpolation = [ux0];
    for i =1:nt_coarse
        %      i = i+1
        U(i+1,:) = BackwardEulerCoarse_heat(dX,dT,f,U(i,:),tC(i),nx_coarse);
        U_interpolation(i+1,:) = spline(xC,U(i+1,:),xF);
    end
    
    % figure
    % surf(U)
    % figure
    % surf(U_interpolation)
    
    % Sequential fine solution
    U_FineSequential=[ux0];
    for n = 1:nt_coarse
        U_Fine = BackwardEulerFine_heat(dx,dt,f,U_FineSequential(n,:),tC(n),nx_fine,m) ;
        U_FineSequential(n+1,:) = U_Fine(end,:);
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
    normInf2executedParearealIter = [];
    normInfDefect = [norm(U_interpolation(end,:) - U_FineSequential(end,:),inf)];
    U0 = U_interpolation;
    
    cmap = hsv(K); % color map
    
    % SupperlinearErrorBound  = normInfDefect;
    % LinearErrorBound  = normInfDefect;
    % proNj = 1;
    
%     figure
%     hold on
%     plot(xF,U_interpolation(end,:),'-*m',xF,U_exact,'-*k');
%     legend('U0','U_{exact}')
%     title(['Solution at time T = ',num2str(T)])
%     xlabel('x')
%     ylabel('u')
    % iteration loop
     alpha = (1-dT);%1;%
        for k = 1: K-1
            % k = k + 1
            % parallel loop - predictor
            for n = 1 : nt_coarse
                % n = n + 1
                F0 =   BackwardEulerFine_heat(dx,dt,f,U_interpolation(n,:),tC(n),nx_fine,m) ;
                ux0F_restriction = restriction(U_interpolation(n,:),nx_coarse,nx_fine);
                G0 =  BackwardEulerCoarse_heat(dX,dT,f,ux0F_restriction,tC(n),nx_coarse)' ;
                G(n,:) = G0(end,:);
                F(n,:) = F0(end,:);
            end % end parallel loop


            % sequential loop - corrector
            if k < K
                for n = 1: nt_coarse
                    %                F0 =  first_order_upwind_adv(a,U(n,:),u0t(n+1),dt,dx,m,'fine');
                    ux0F_restriction = restriction(U_interpolation(n,:),nx_coarse,nx_fine);
                    Gn = BackwardEulerCoarse_heat(dX,dT,f,ux0F_restriction,tC(n),nx_coarse)' ;
                    U_interpolation( n + 1,:) = alpha*(spline(xC,Gn(end,:),xF)  - spline(xC,G(n,:),xF)) + F(n,:);
                    %                U( n + 1,:) = restriction(U_temp,nx_coarse);
                    norm((spline(xC,Gn(end,:),xF)  - spline(xC,G(n,:),xF)))
                end % sequential loop
            end
            normInf_U0_Uinterpolate = norm(U0(end,:)-U_interpolation(end,:),inf) ;
            normInf2executedParearealIter =[normInf2executedParearealIter;normInf_U0_Uinterpolate] ;
            normInfDefect =[normInfDefect ; norm(U_interpolation(end,:)-U_FineSequential(end,:),inf)];

            U0 = U_interpolation;
    %         plot(xF,U_interpolation(end,:),'Color',cmap(k,:),'Marker','o');

    
        end  % iteration loop
    normInf2executedParearealIter
    normInfDefect

    
    figure
    semilogy(1:K-1, normInf2executedParearealIter,'-*b',0:K-1,normInfDefect,'-*r')
    legend('normInf2executedParearealIter','normInfDefect')
    xlabel('k')
    title('BackwardEuler')
end



