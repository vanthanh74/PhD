function [x,flag,relres,iter,resvec,L2NormErrorGMRES,LInfNormErrorGMRES] = mygmresODE(A,b,restart,tol,maxit,M1,M2,x,varargin)
%GMRES   Generalized Minimum Residual Method.
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

global a N m K T;
a = -1;
yExactFunction = @(t)exp(a*t);

% stability function (  y^n+1 = R(z)y^n  )
Rz = @(z) 1/(1-z);     % Backward Euler
% Rz = @(z) 1+z;       % Forward Euler

%% -----------------------------------------------------------------


% initial condition
u0 = 1;

% coarse propagator
% N = 40;     %   number of coarse time intervals
% T = 1;      %   T = N*Dt
Dt = T/N;   %   coarse timesteps

xt = linspace(0,T,N+1);

% fine propagator
% m = 10;    % number of fine time steps in each coarse time interval
dt = Dt/m;

% K = 40; % number of Parareal iterations

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


% coarse space initial
indice_R_coarse  = 0:m:m*N;
R = zeros(length(indice_R_coarse),N_fine_nodes);
for i = 1:length(indice_R_coarse)
    R(i,indice_R_coarse(i)+1) = 1;
end

% first subdomain  Omega_1 = {0}
R1 = zeros(1,N_fine_nodes);
R1(1) = 1;
A1 = R1*A_original*R1';

% define nodes in other subdomains
N_nodes_subdomain = m;
indice_R_subdomains = zeros(N,N_nodes_subdomain);
for i =1:N
    indice_R_subdomains(i,:) =  (i-1)*m+1:(i-1)*m+m;
end

% Omega_2 = {1,2} Omega_3 = {3,4}
R2_end = zeros(N*N_nodes_subdomain,N_fine_nodes) ;
for i = 1:N
    for j = 1:N_nodes_subdomain
        R2_end(N_nodes_subdomain*(i-1)+j,indice_R_subdomains(i)+j) = 1;
    end
end

R_sub = [];       % restriction matrices
A_sub = [];       % sub-matrices
for i = 2: N+1
    R_sub{i} = R2_end(N_nodes_subdomain*(i-2)+1:N_nodes_subdomain*(i-2)+N_nodes_subdomain,:);
    A_sub{i} = R_sub{i}*A_original*R_sub{i}';
end

% define 2-level domain decomposition multiplicative preconditioner
P = R1'*A1*R1;
for i = 2:N+1
    P = P + R_sub{i}'*A_sub{i}*R_sub{i};
end

%Multiplicative
M = P*(R'*B*R +  Id-R'*R);


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

L2NormErrorGMRES = sqrt(sum(Dt*(uF_true(pick)-x(pick)).^2));%norm(uF_true(pick)-U_fine(pick),2);%
LInfNormErrorGMRES = norm(uF_true(pick)-x(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));
SupperlinearErrorBound  = L2NormErrorGMRES;
LinearErrorBound  = LInfNormErrorGMRES;

% % rhs
% f = zeros(N_fine_nodes,1);
% f(1) = u0;




















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
%         maxit = min(n,10);
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
resvec(1) =  1;% relative residual ||r||/||r0||%norm(f-A_original*U0_fine);% normr;%               % resvec(1) = norm(b-A*x0)
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
    
%         proNj(initer+1) = proNj(initer)*(N-initer);
%         L2NormErrorGMRES(initer+1) = sqrt(sum(Dt*(uF_true(pick)-XM(pick)).^2));%norm(uF_true(pick)-U_fine(pick),2);%
%         LInfNormErrorGMRES(initer+1) = norm(uF_true(pick)-XM(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));%
%         SupperlinearErrorBound(initer+1)  = (abs(exp(a*Dt) - Rz(a*Dt)))^(initer)/factorial(initer)*proNj(initer+1)*max(abs(uF_true(pick)-U0n(pick)));
%         LinearErrorBound(initer+1)  = (abs(exp(a*Dt) - Rz(a*Dt))/(1- abs(Rz(a*Dt))))^(initer)*max(abs(uF_true(pick)-U0n(pick)));
        
        
        
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
        
        
        proNj(initer+1) = proNj(initer)*(N-initer);
        L2NormErrorGMRES(initer+1) = sqrt(sum(Dt*(uF_true(pick)-XM(pick)).^2));%norm(uF_true(pick)-U_fine(pick),2);%
        LInfNormErrorGMRES(initer+1) = norm(uF_true(pick)-XM(pick),inf);%max(abs(uF_true(pick)-U_fine(pick)));%
        SupperlinearErrorBound(initer+1)  = (abs(exp(a*Dt) - Rz(a*Dt)))^(initer)/factorial(initer)*proNj(initer+1)*max(abs(uF_true(pick)-U0n(pick)));
        LinearErrorBound(initer+1)  = (abs(exp(a*Dt) - Rz(a*Dt))/(1- abs(Rz(a*Dt))))^(initer)*max(abs(uF_true(pick)-U0n(pick)));
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
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

