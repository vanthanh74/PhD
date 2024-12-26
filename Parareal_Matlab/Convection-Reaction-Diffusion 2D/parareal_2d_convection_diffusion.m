function  [nmv_total_matrix_FineSolver,nmv_total_matrix_FineSolverInParareal,TimeConsumingFineSolver,TimeConsumingParareal,total_nmv_perIterParareal,total_nmv] = parareal_2d_convection_diffusion(L,T,kappa,dT,dX_coarse,dY_coarse,...
    dt,dx_fine,dy_fine,nt_coarse,nx_coarse,ny_coarse,...
    nt_fine,nx_fine,ny_fine,m,K)


% IVP: u_t = kappa*(u_xx + u_yy) - 2*p1*u_x - 2*p2*u_y + p3*u + f    in [0,L]x[0,L] x (0,T),
%      u(0,y,t) = u(L,y,t)  = 0                 0 <= y <= L, t >= 0
%      u(x,0,t) = u(x,L,t)  = 0                 0 <= x <= L, t >= 0
%      u(x,y,0) = u0(x,y)                   (x,y) in [0,L]x[0,L]

% INPUT PARAMETERS:
%   * nt_coarse: number of coarse intervals in time
%   * nx_coarse, ny_coarse: number of coarse intervals in space
%   * nt_fine: number of fine intervals in time
%   * nx_fine, ny_fine: number of fine intervals in space
%   * m: number of fine time steps on each coarse time step 
%   * K: number of Parareal iterations (2 <= K <= n_coarse+1)
%        (--> K sequential coarse sweeps and K-1 parallel fine sweeps)
% OUTPUT:
%   * LInfinityErrorNorm


%% initial data
global xC yC xF yF p1 p2 p3

% coarse points
xC = linspace(0,L,nx_coarse+1);
yC = linspace(0,L,ny_coarse+1);
tC = linspace(0,T,nt_coarse+1);

% fine points
xF = linspace(0,L,nx_fine+1);
yF = linspace(0,L,ny_fine+1);
tF = linspace(0,T,nt_fine+1);

% exact solution
u_ex=@(x,y,t) x.*(L-x).^2.*y.*(L-y).^2.*exp(-2*t);

% Initial condition u(x,y,0) = u0(x)
u0_function = @(x,y) x.*(L-x).^2.*y.*(L-y).^2;
% u0_function = @(x,y) sin(2*pi*x);
% u0_function = @(x,y) 0;




% f = diff(u,t) - kappa*(diff(diff(u,x),x) + diff(diff(u,y),y)) +2*p1*diff(u,x) +2*p2*diff(u,y) -p3*u 

% function f

 
% f = @(x,y,t) -2* exp(-2* t) *(L - y)* y* (L - x)* x...
%          - (-2* exp(-2 *t) *(L - y)* y - 2* exp(-2* t)* (L - x)* x) *kappa...
%          + 2 *(-exp(-2 *t)* (L - y)* y *x + exp(-2* t)* (L - y)* y* (L - x))* p1...
%          + 2* (-exp(-2* t)* y *(L - x)* x + exp(-2* t)* (L - y) *(L - x)* x) *p2...
%          - exp(-2* t)* (L - y) *y *(L - x) *x* p3;

         
                  
  f = @(x,y,t)  -2*exp(-2*t)*(L - y)^2*y*(L - x)^2* x - (2*exp(-2*t)*(L - y)^2*y*x ...                                                              
                - 4*exp(-2*t)*(L - y)^2*y*(L - x) + 2*exp(-2*t)*y*(L - x)^2*x  ...                                   
                  - 4*exp(-2*t)*(L - y)*(L - x)^2*x)*kappa  ...
                  + 2*(-2*exp(-2*t)*(L - y)^2*y*(L - x)*x + exp(-2*t)*(L - y)^2*y*(L - x)^2 )*p1...
                 + 2*(-2*exp(-2*t)*(L - y)*y*(L - x)^2*x + exp(-2*t)*(L - y)^2*(L - x)^2*x)*p2   ...                     
                  - exp(-2*t)*(L - y)^2*y*(L - x)^2*x*p3;

% f = @(x,y,t) sqrt(x*y*t)*sin(2*x*y*t);

% f = @(x,y,t) x^4*(1-x)+t^2;

% initial solution u0
% u0 = zeros(1,(nx_fine-1)*(ny_fine-1));
% for i=1:nx_fine-1
%     for j=1:ny_fine-1
%         u0((i-1)*(nx_fine-1) +j)= u0_function(i*dx_fine,j*dy_fine);
%     end
% end

u0 = zeros(nx_fine+1,ny_fine+1);
for i=1:nx_fine+1
    for j=1:ny_fine+1
        u0(i,j)= u0_function(xF(i),yF(j));
    end
end
%% ------------------------------------------------------------------
% initial fine integration
% restriction
u0_restriction = restriction(u0,nx_coarse,ny_coarse,nx_fine,ny_fine);
%     U = [u0_restriction];
U_interpolation = [u0];

for i =1 :nt_coarse
    U = BackwardEulerCoarse_2DConvDiff(kappa,dX_coarse,dY_coarse,dT,f,u0_restriction,tC(i),nx_coarse,ny_coarse);
    u0_restriction = U;
    
    U_interp_temp = interpolation(xC,yC,U,xF,yF);
    %     figure
    %     surf(xF,yF,U_interp_temp)
    %          U_interpolation((i-1)*(nx_fine-1)+1:(i)*(nx_fine-1),(i-1)*(ny_fine-1)+1:(i)*(ny_fine-1))=U_interp_temp;
    U_interpolation = [ U_interpolation ; U_interp_temp];
    
end




% Sequential fine solution
name = strseq('subspace',1:nt_coarse);
name = char(name);
U_FineSequential=[u0];
U_Fine0 = u0;
nmv_total_matrix_FineSolver = []; %Record number of matrix-vector products needed to solve this linear system
tic
for n = 1:nt_coarse
    if nt_coarse >= 10 && nt_coarse <= 99 && n <=9
           [U_Fine,nmv_total] = BackwardEulerFine_2DConvDiff(kappa,dx_fine,dy_fine,dt,f,U_Fine0,tC(n),nx_fine,ny_fine,m,name(n,1:end-1));
    elseif nt_coarse >= 10 && nt_coarse <= 99 && n >=10
           [U_Fine,nmv_total] = BackwardEulerFine_2DConvDiff(kappa,dx_fine,dy_fine,dt,f,U_Fine0,tC(n),nx_fine,ny_fine,m,name(n,:));
    end
    U_Fine0 = U_Fine(end-nx_fine:end,:);
    U_FineSequential = [U_FineSequential;U_Fine(end-nx_fine:end,:) ];
    nmv_total_matrix_FineSolver = [nmv_total_matrix_FineSolver, nmv_total];
end
TimeConsumingFineSolver = toc;


%  nmv_total_matrix_FineSolver
%   nmv_total_matrix_FineSolverInParareal = [];
%  TimeConsumingParareal = [];
%  return
 
 
% Exact solution at time T = 1
U_exact=zeros(nx_fine+1,ny_fine+1);
for i=1:nx_fine+1
    for j=1:ny_fine+1
        U_exact(i,j)=u_ex(xF(i),yF(j),T);
    end
end

norm(U_FineSequential(end-nx_fine:end,:) - U_exact,inf);
norm(U_interpolation(end-nx_fine:end,:)  - U_exact,inf);

normInf2executiveParearealIter = [];
normInfDefect = [norm(U_interpolation(end-nx_fine:end,:) - U_FineSequential(end-nx_fine:end,:),inf)];
U0 = U_interpolation;


figure
subplot(121)
surf(xF,yF,U_interpolation(end-nx_fine:end,:))
zlim([0 3*1e-3])
xlabel('x')
ylabel('y')
zlabel('u')
title('U0')

subplot(122)
surf(xF,yF,U_exact)
zlim([0 3*1e-3])
xlabel('x')
ylabel('y')
zlabel('u')
title('U_{exact}')



cmap = hsv(K); % color map
nmv_total_matrix_FineSolverInParareal =[];
nameFineSolver = strseq('subspaceFine',1:nt_coarse*(K-1));
nameFineSolver = char(nameFineSolver);
total_nmv_perIterParareal =[] ; % number of matrix vector products per each Parareal iteration
tic
% iteration loop
for k = 1: K-1
    k
    % k = k + 1
    G = [];
    F = [];
    
    % parallel loop - predictor
    for  n = 1 : nt_coarse
%          n = n + 1
%          nameFineSolver(nt_coarse*(k-1) + n,:)
%          nameFineSolver(nt_coarse*(k-1) + n,1:end-1)
%         nameFineSolver(nt_coarse*(k-1) + n,1:end-2)
        if nt_coarse*(k-1) + n <= 9
            [F0 ,nmv_toltal]=   BackwardEulerFine_2DConvDiff(kappa,dx_fine,dy_fine,dt,f,U_interpolation((n-1)*(nx_fine+1)+1:n*(nx_fine+1),:),tC(n),nx_fine,ny_fine,m,nameFineSolver(nt_coarse*(k-1) + n,1:end-2));
        elseif nt_coarse*(k-1) + n >=10 && nt_coarse*(k-1) + n <= 99
            [F0 ,nmv_toltal]=   BackwardEulerFine_2DConvDiff(kappa,dx_fine,dy_fine,dt,f,U_interpolation((n-1)*(nx_fine+1)+1:n*(nx_fine+1),:),tC(n),nx_fine,ny_fine,m,nameFineSolver(nt_coarse*(k-1) + n,1:end-1));
        elseif nt_coarse*(k-1) + n >=100 && nt_coarse*(k-1) + n <= 999
            [F0 ,nmv_toltal]=   BackwardEulerFine_2DConvDiff(kappa,dx_fine,dy_fine,dt,f,U_interpolation((n-1)*(nx_fine+1)+1:n*(nx_fine+1),:),tC(n),nx_fine,ny_fine,m,nameFineSolver(nt_coarse*(k-1) + n,:));            
         end
            ux0F_restriction = restriction(U_interpolation((n-1)*(nx_fine+1)+1:n*(nx_fine+1),:),nx_coarse,ny_coarse,nx_fine,ny_fine);
            G0 =   BackwardEulerCoarse_2DConvDiff(kappa,dX_coarse,dY_coarse,dT,f,ux0F_restriction,tC(n),nx_coarse,ny_coarse) ;
        G = [G;G0];
        F = [F;F0(end-nx_fine:end,:)];    
         total_nmv_perIterParareal(k,n)=  sum(nmv_toltal);
        if k==1
            nmv_total_matrix_FineSolverInParareal((k-1)*m+1:k*m,n) = nmv_toltal;
        else
            nmv_total_matrix_FineSolverInParareal((k-1)*m+k:k*m+k-1,n) = nmv_toltal;
        end
        if n == nt_coarse
            nmv_total_matrix_FineSolverInParareal(end+1,:) = zeros(1,nt_coarse);
        end
    end
    %     size(G)
    %     size(F)
   
        
%     TimeConsumingParareal=[]
%     return
       
     % sequential loop - corrector
    if k < K
        for n = 1: nt_coarse
            ux0F_restriction = restriction(U_interpolation((n-1)*(nx_fine+1)+1:n*(nx_fine+1),:),nx_coarse,ny_coarse,nx_fine,ny_fine);
            Gn = BackwardEulerCoarse_2DConvDiff(kappa,dX_coarse,dY_coarse,dT,f,ux0F_restriction,tC(n),nx_coarse,ny_coarse) ;
            U_interpolation( (n)*(nx_fine+1)+1:(n+1)*(nx_fine+1),:) =     interpolation(xC,yC,Gn,xF,yF)  ...
                - interpolation(xC,yC,G((n-1)*(nx_coarse+1)+1:n*(nx_coarse+1),:),xF,yF) ...
                + F((n-1)*(nx_fine+1)+1:n*(nx_fine+1),:);
        end
    end
    
    
    normInf_U0_Uinterpolate = norm(U0(end-nx_fine:end,:)-U_interpolation(end-nx_fine:end,:),inf) ;
    normInf2executiveParearealIter =[normInf2executiveParearealIter;normInf_U0_Uinterpolate] ;
    normInfDefect =[normInfDefect ; norm(U_interpolation(end-nx_fine:end,:)-U_FineSequential(end-nx_fine:end,:),inf)];
    U0 = U_interpolation;
end
normInf2executiveParearealIter
normInfDefect
TimeConsumingParareal=toc;


figure
surf(xF,yF,real(U_interpolation(end-nx_fine:end,:)))
zlim([0 3*1e-3])
xlabel('x')
ylabel('y')
zlabel('u')
title(['Parareal solution at T = ',num2str(T)])

figure
surf(xF,yF,real(U_FineSequential(end-nx_fine:end,:)))
zlim([0 3*1e-3])
xlabel('x')
ylabel('y')
zlabel('u')
title(['Sequential solution at T = ',num2str(T)])

figure
surf(xF,yF,U_exact)
zlim([0 3*1e-3])
xlabel('x')
ylabel('y')
zlabel('u')
title(['Exact solution at T = ',num2str(T)])

figure
subplot(131)
surf(xF,yF,real(U_FineSequential(end-nx_fine:end,:)))
zlim([0 3*1e-3])
xlabel('x')
ylabel('y')
zlabel('u')
title('U_{fine}')

subplot(132)
surf(xF,yF,U_exact)
zlim([0 3*1e-3])
xlabel('x')
ylabel('y')
zlabel('u')
title('U_{exact}')

subplot(133)
surf(xF,yF,real(U_interpolation(end-nx_fine:end,:)))
zlim([0 3*1e-3])
xlabel('x')
ylabel('y')
zlabel('u')
title('Parareal solution')

figure
semilogy(1:K-1, normInf2executiveParearealIter,'-*b',0:K-1,normInfDefect,'-*r')
legend('normInf2executiveParearealIter','normInfDefect')
xlabel('k')
title('BackwardEuler')

total_nmv = sum(sum(total_nmv_perIterParareal));

