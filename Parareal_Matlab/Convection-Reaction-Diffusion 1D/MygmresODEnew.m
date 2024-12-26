function [x, error, res, iter, flag, L2NormErrorGMRES, LInfNormErrorGMRES] = MygmresODEnew( A, x, b, M, restrt, max_it, tol )

%  -- Iterative template routine --
%     Univ. of Tennessee and Oak Ridge National Laboratory
%     October 1, 1993
%     Details of this algorithm are described in "Templates for the
%     Solution of Linear Systems: Building Blocks for Iterative
%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
%     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
%
% [x, error, iter, flag] = gmres( A, x, b, M, restrt, max_it, tol )
%
% gmres.m solves the linear system Ax=b
% using the Generalized Minimal residual ( GMRESm ) method with restarts .
%
% input   A        REAL nonsymmetric positive definite matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         M        REAL preconditioner matrix
%         restrt   INTEGER number of iterations between restarts
%         max_it   INTEGER maximum number of iterations
%         tol      REAL error tolerance
%
% output  x        REAL solution vector
%         error    REAL error norm
%         iter     INTEGER number of iterations performed
%         flag     INTEGER: 0 = solution found to tolerance
%                           1 = no convergence given max_it




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
N = 20;     %   number of coarse time intervals
T = 1;      %   T = N*Dt
Dt = T/N;   %   coarse timesteps

xt = linspace(0,T,N+1);

% fine propagator
m = 2;    % number of fine time steps in each coarse time interval
dt = Dt/m;

K = 20; % number of Parareal iterations

cmap = hsv(K); % color map

% 2-level domain decomposition initial
N_fine_nodes = m*N+1; % 11
% N_subdomains = m+1; % 6
Id = eye(N_fine_nodes);

% plot the time domain

% figure
% plot(0:N_fine_nodes-1,1,'.b',0:m:N_fine_nodes-1,1,'ro')
% xlim([0 N_fine_nodes-1])
% ylim([0.5 1.5])

%%
% original problem matrix
phi = (1-a*dt)^-1;
A_original = eye(m*N+1,m*N+1) + (-phi)*diag(ones(m*N,1),-1);

% Reduced system

% fine propagator
F_tilde = (1-a*dt)^-m;

% coarse propagator
G_tilde = (1-a*Dt)^-1;



% Uk_0 = zeros(m*N,1);
% Uk_0(1) = u0;

rhs0 = [u0; zeros(m*N,1)];
% rhs = [u0; phi*u0;-G_tilde*u0;phi*U(3);-G_tilde*U(3)];


% exact solution

ttrue = linspace(0,T,N*m+1);
pick = 1:m:N*m+1;

uF_true = A_original\rhs0;

% reduced matrix
A_tilde = eye(N+1,N+1) + (-F_tilde)*diag(ones(N,1),-1);
% A_tilde = eye(m*N+1,m*N+1) + (-F_tilde)*diag(ones(m*N,1),-1);

% approximate reduced matrix
B = eye(N+1,N+1) + (-G_tilde)*diag(ones(N,1),-1);

% initial coarse integration

B1 =  eye(N+1,N+1) + (-G_tilde)*diag(ones(N,1),-1);
% B1 =  eye(N*m+1,N*m+1) + (-G_tilde)*diag(ones(N*m,1),-1);

%Initial coarse Parareal iteration at the coarse level
rhs0_coarse = [u0;zeros(N,1)];
U0_coarse = B1\rhs0_coarse;

U_fine=spline(ttrue(pick),U0_coarse,ttrue)';
U0_fine = U_fine;

% figure
% hold on
% plot(ttrue,uF_true,'-^k',ttrue,U_fine,'-*m');
% legend('True solution','Parareal solution at k = 0')

U0n = U_fine;
proNj = 1;

L2NormErrorGMRES = sqrt(sum(Dt*(uF_true(pick)-U_fine(pick)).^2));%norm(uF_true(pick)-U_fine(pick),2);%
LInfNormErrorGMRES = norm(uF_true(pick)-U_fine(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));
SupperlinearErrorBound  = L2NormErrorGMRES;
LinearErrorBound  = LInfNormErrorGMRES;

% rhs
f = zeros(N_fine_nodes,1);
f(1) = u0;


%%
%------------------------------------------------------------------------------------
   iter = 0;                                         % initialization
   flag = 0;
   x0 = x;
   bnrm2 = norm( b );
   if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end

   r = M \ ( b-A*x );
   r0 = r;
   res = norm(r)/norm(r0);
   error = norm( r )/bnrm2 ;
   
   if ( error < tol ) return, end

   [n,n] = size(A);                                  % initialize workspace
   m = restrt;
   V(1:n,1:m+1) = zeros(n,m+1);
   H(1:m+1,1:m) = zeros(m+1,m);
   cs(1:m) = zeros(m,1);
   sn(1:m) = zeros(m,1);
   e1    = zeros(n,1);
   e1(1) = 1.0;

   for iter = 1:max_it                              % begin iteration

      r = M \ ( b-A*x );
      if iter > 1
           res = norm(r)/norm(r0);
      end
      V(:,1) = r / norm( r );
      s = norm( r )*e1;
      for i = 1:m                                   % construct orthonormal
	 w = M \ (A*V(:,i));                         % basis using Gram-Schmidt
	 for k = 1:i
	   H(k,i)= w'*V(:,k);
	   w = w - H(k,i)*V(:,k);
	 end
	 H(i+1,i) = norm( w );
	 V(:,i+1) = w / H(i+1,i);
	 for k = 1:i-1                              % apply Givens rotation
            temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i);
            H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
            H(k,i)   = temp;
	 end
	 [cs(i),sn(i)] = rotmat( H(i,i), H(i+1,i) ); % form i-th rotation matrix
         temp   = cs(i)*s(i);                        % approximate residual norm
         s(i+1) = -sn(i)*s(i);
	 s(i)   = temp;
         H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
         H(i+1,i) = 0.0;
	 error  = abs(s(i+1)) ;
%      proNj(iter+1) = proNj(iter)*(N-iter);
%      L2NormErrorGMRES(iter+1) = sqrt(sum(Dt*(uF_true(pick)-x(pick)).^2));%norm(uF_true(pick)-U_fine(pick),2);%
%      LInfNormErrorGMRES(iter+1) = norm(uF_true(pick)-x(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));%
%      SupperlinearErrorBound(iter+1)  = (abs(exp(a*Dt) - Rz(a*Dt)))^(iter)/factorial(iter)*proNj(iter+1)*max(abs(uF_true(pick)-U0n(pick)));
%      LinearErrorBound(iter+1)  = (abs(exp(a*Dt) - Rz(a*Dt))/(1- abs(Rz(a*Dt))))^(iter)*max(abs(uF_true(pick)-U0n(pick)));

        
	 if ( error <= tol )                        % update approximation
         y = H(1:i,1:i) \ s(1:i);                 % and exit
         x = x + V(:,1:i)*y
         r = M \ ( b-A*x )                              % compute residual
         res(iter+1,:) = norm(r)/norm(r0);
         s(i+1) = norm(r);
         error = s(i+1) ;         % check convergence
         % Compute error norms
         proNj(iter+1) = proNj(iter)*(N-iter);
         L2NormErrorGMRES(iter+1) = sqrt(sum(Dt*(uF_true(pick)-x(pick)).^2));%norm(uF_true(pick)-U_fine(pick),2);%
         LInfNormErrorGMRES(iter+1) = norm(uF_true(pick)-x(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));%
         SupperlinearErrorBound(iter+1)  = (abs(exp(a*Dt) - Rz(a*Dt)))^(iter)/factorial(iter)*proNj(iter+1)*max(abs(uF_true(pick)-U0n(pick)));
         LinearErrorBound(iter+1)  = (abs(exp(a*Dt) - Rz(a*Dt))/(1- abs(Rz(a*Dt))))^(iter)*max(abs(uF_true(pick)-U0n(pick)));

        	    break;
	 end
      end

        
      if ( error <= tol ), break, end
      y = H(1:m,1:m) \ s(1:m);
      x = x + V(:,1:m)*y                           % update approximation
      r = M \ ( b-A*x )                              % compute residual
      res(iter+1,:) = norm(r)/norm(r0);
      s(i+1) = norm(r);
      error = s(i+1) ;         % check convergence
             proNj(iter+1) = proNj(iter)*(N-iter);
         L2NormErrorGMRES(iter+1) = sqrt(sum(Dt*(uF_true(pick)-x(pick)).^2));%norm(uF_true(pick)-U_fine(pick),2);%
         LInfNormErrorGMRES(iter+1) = norm(uF_true(pick)-x(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));%
         SupperlinearErrorBound(iter+1)  = (abs(exp(a*Dt) - Rz(a*Dt)))^(iter)/factorial(iter)*proNj(iter+1)*max(abs(uF_true(pick)-U0n(pick)));
         LinearErrorBound(iter+1)  = (abs(exp(a*Dt) - Rz(a*Dt))/(1- abs(Rz(a*Dt))))^(iter)*max(abs(uF_true(pick)-U0n(pick)));

      if ( error <= tol ), break, end
      
%        proNj(iter+1) = proNj(iter)*(N-iter);
%          L2NormErrorGMRES(iter+1) = sqrt(sum(Dt*(uF_true(pick)-x(pick)).^2));%norm(uF_true(pick)-U_fine(pick),2);%
%          LInfNormErrorGMRES(iter+1) = norm(uF_true(pick)-x(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));%
%          SupperlinearErrorBound(iter+1)  = (abs(exp(a*Dt) - Rz(a*Dt)))^(iter)/factorial(iter)*proNj(iter+1)*max(abs(uF_true(pick)-U0n(pick)));
%          LinearErrorBound(iter+1)  = (abs(exp(a*Dt) - Rz(a*Dt))/(1- abs(Rz(a*Dt))))^(iter)*max(abs(uF_true(pick)-U0n(pick)));
	   
  
    
      
      
      
   end

   if ( error > tol ) flag = 1; end                 % converged
%    res = s;
% END of gmres.m