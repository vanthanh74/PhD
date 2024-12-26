function DU = UpdateConvDiff(FF,U_old)

global xC xF T L a bb c nt_interval nx_coarse nt_coarse nx_fine nt_fine

% T = 1; %  Intervall (0,T)
% L = 1; % omega=(0,L)

%--------------------------------------
% % % solution x(L-x)^2*exp(-2*t);

% a = 3;%3
% b = 0;%0.005
% c = 0;%1

% nt_interval = 20; % number of time intervals
% 
% nx_coarse = 10; % number of coarse spamtial intervals
% nt_coarse = 20;% 300 %20    % number of coarse time intervals
% 
% 
% nx_fine = 	10; % number of fine spatial intervals
% nt_fine = 40;% 3000 %200 . % number of fine time intervals (nt_coarsex2 )

%%%%--------------------------------------
% figure
% spy(P)
% figure
% spy(R'*M_tilde*R + Id - R'*R)
% figure
% spy(M_mul)

%--------------------------------------
% solution sin(2*pi*x)*exp(-2*t);

% a = 0.5;      %0.5;
% b = 0;      %0.0025;
% c = 0;      %0;
% 
% nt_interval = 20; % number of time intervals
% 
% nx_coarse = 10;
% nt_coarse = 20;% 300 %20    % number of coarse time intervals
% 
% 
% nx_fine = 	10;
% nt_fine = 40;% 3000 %200 . % number of fine time intervals



%--------------------------------------
% % % % solution x*(L-x)*exp(-2*t);
%
% a = 3;
% b = 0.005;
% c = 1;
%
% nt_interval = 10; % number of time intervals
%
% nx_coarse = 10; % number of coarse spatial intervals
% nt_coarse = 10;% 300 %20    % number of coarse time intervals
%
%
% nx_fine = 	10; % number of fine spatial intervals
% nt_fine = 160;% 3000 %200 . % number of fine time intervals (nt_coarsex2 )


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

M = DT/dT ; % number of coarse time steps on each time interval
m = DT/dt ; % number of fine time steps on each time interval

r = dT/dX^2;
r1 = dt/dx^2;

K = 20; % number of Parareal iterations

% 2-level domain decomposition initial
N_fine_nodes = m*nt_interval+1 ;% 11
% N_subdomains = m + 1% nt_interval + 1 % 6
% Id = eye(N_fine_nodes);


% % plot the time domain
% figure
% plot(0:N_fine_nodes-1,1,'.b',0:m:N_fine_nodes-1,1,'ro')
% xlim([0 N_fine_nodes-1])
% ylim([0.5 1.5])

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


%-------------------------------------------------------
% % exact solution 1
u_ex=@(x,t) x*(L-x)^2*exp(-2*t);%sin(2*pi*x)*exp(-2*t);%

ux0_function = @(x) (x.*(L-x).^2);%sin(2*pi*x);%

f = @(x,t) bb.*(exp(-2.*t).*(L - x).^2 - x.*exp(-2.*t).*(2.*L - 2.*x)) - a.*(2.*x.*exp(-2.*t) - 2.*exp(-2.*t).*(2.*L - 2.*x)) - 2.*x.*exp(-2.*t).*(L - x).^2 - c.*x.*exp(-2.*t).*(L - x).^2;
% % f = @(x,t) x.^4.*(1-x)+t.^2;


%-------------------------------------------------------
% % exact solution 2
% u_ex=@(x,t) sin(2*pi*x)*exp(-2*t);%x*(L-x)^2*exp(-2*t);
% 
% %syms a b c x t L
% % u = sin(2*pi*x)*exp(-2*t);%x*(L-x)^2*exp(-2*t)
% % f = diff(u,t) - a*diff(diff(u,x),x) + b*diff(u,x) - c*u
% 
% ux0_function = @(x) sin(2*pi*x);%(x.*(L-x).^2);
% 
% f = @(x,t) 4*a*pi^2*exp(-2*t)*sin(2*pi*x) - c*exp(-2*t)*sin(2*pi*x) - 2*exp(-2*t)*sin(2*pi*x) + 2*b*pi*exp(-2*t)*cos(2*pi*x);

%-------------------------------------------------------
% % exact solution 3
% % % syms a b c x t L
% % % u = x*(L-x)*exp(-2*t)
% % % f = diff(u,t) - a*diff(diff(u,x),x) + b*diff(u,x) - c*u
% u_ex=@(x,t) x*(L-x)*exp(-2*t);%sin(2*pi*x)*exp(-2*t);%
%
% ux0_function = @(x) (x.*(L-x));%sin(2*pi*x);%
%
% f = @(x,t) 2*a*exp(-2*t) - b*(x.*exp(-2*t) - exp(-2*t).*(L - x)) - 2*x.*exp(-2*t).*(L - x) - c*x.*exp(-2*t).*(L - x);



%% Define fine and coarse propagators
% Backward Euler fine solver
% A_tilde
rF = a*dt/(dx^2);
kF = bb*dt/(2*dx);


A_F = (kF-rF)*diag(ones(nx_coarse-2,1),1) + (1+2*rF-c*dt)*eye(nx_coarse-1) - (kF+rF)*diag(ones(nx_coarse-2,1),-1);
Id_F = eye(size(A_F));
% A_F_Global =  zeros((nt_interval+1)*(nx_coarse-1));
% for i =1:nt_interval+1
%     %  i = i + 1
%     A_F_Global((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i-1)*(nx_coarse-1)+1:i*(nx_coarse-1)) = Id_F;
%     if i>1
%         A_F_Global((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i-2)*(nx_coarse-1)+1:(i-1)*(nx_coarse-1)) = -A_F^-(m);
%     end
% end
% A_tilde = A_F_Global;
% condNumberA_F = cond(-A_F^-(m));
% condNumberA_F_Global=cond(A_F_Global);
% Backward Euler coarse solver
% M_tilde
rC = a*dT/(dX^2);
kC = bb*dT/(2*dX);

A_C = (kC-rC)*diag(ones(nx_coarse-2,1),1) + (1+2*rC-c*dT)*eye(nx_coarse-1) - (kC+rC)*diag(ones(nx_coarse-2,1),-1);
Id_C = eye(size(A_C));
A_C_Global =  zeros((nt_interval+1)*(nx_coarse-1));
for i =1:nt_interval+1
    %  i = i + 1
    A_C_Global((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i-1)*(nx_coarse-1)+1:i*(nx_coarse-1)) = Id_C;
    if i>1
        A_C_Global((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i-2)*(nx_coarse-1)+1:(i-1)*(nx_coarse-1)) = -A_C^-1;
    end
end
M_tilde = A_C_Global;
% condNumberA_C = cond(-A_C^-(1))
% condNumberA_C_Global = cond(A_C_Global)
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
% U_exact=zeros(1,nx_fine+1);
% U_exact_full = zeros((nt_interval+1)*(nx_fine-1),1);
% for i=1:nx_fine+1
%     U_exact(i)=u_ex(xF(i),T);
% end

% initial solution u0
ux0 = ux0_function(xF(2:end-1))';

% Fine sequential propagator
% Uk_0 = zeros((nt_interval+1)*(nx_fine-1),1);
% Uk_0(1:nx_fine-1) = ux0;
% rhs0 = [ux0; zeros(((nt_interval+1)*(nx_fine-1)-nx_fine+1),1)];
%  rhs = [ux0; (A_tilde-M_tilde)*Uk_0];

% % Define original A and I_F
A_Fo = (kF-rF)*diag(ones(nx_fine-2,1),1) + (1+2*rF-c*dt)*eye(nx_fine-1) - (kF+rF)*diag(ones(nx_fine-2,1),-1);
Id_Fo = eye(size(A_Fo));
A_Original =  zeros(m*nt_interval+1);
for i =1:m*nt_interval+1
    %  i = i + 1
    A_Original((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i-1)*(nx_coarse-1)+1:i*(nx_coarse-1)) = Id_Fo;
    if i>1
        A_Original((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i-2)*(nx_coarse-1)+1:(i-1)*(nx_coarse-1)) = -A_Fo^-1;
    end
end
% condNumber_A_Original=cond(A_Original)



% % Restriction matric for the right hand side
% R_rhs =  zeros((nt_interval+1)*(nx_fine-1),(m*nt_interval+1)*(nx_fine-1));
% seq_rhs = [zeros(nx_fine-1,(m*nt_interval+1)*(nx_fine-1))];
% for i = 1:m
%     seq_rhs(:,(i)*(nx_fine-1)+1:(i+1)*(nx_fine-1)) =  A_Fo^-(m-i);
% end
% 
% for i = 1:nt_interval+1
%     %  i = i + 1
%     if i == 1
%         R_rhs((i-1)*(nx_fine-1)+1:i*(nx_fine-1),(i-1)*(nx_fine-1)+1:i*(nx_fine-1)) = Id_Fo;
%         i_temp = 0;
%         
%     else
%         %         R_rhs((i-1)*(nx_fine-1)+1:i*(nx_fine-1),(i-1+i_temp)*(nx_fine-1)+1:(i+i_temp)*(nx_fine-1)) = A_F^-1;
%         %         R_rhs((i-1)*(nx_fine-1)+1:i*(nx_fine-1),(i+i_temp)*(nx_fine-1)+1:(i+1+i_temp)*(nx_fine-1)) = Id_F;
%         %         i_temp = i_temp + 1;
%         
%         R_rhs((i-1 )*(nx_fine-1)+1:(i )*(nx_fine-1),:) = seq_rhs;
%         seq_rhs = [zeros(nx_fine-1,(m)*(nx_fine-1)) seq_rhs];
%         seq_rhs(:,end-(m)*(nx_fine-1)+1:end)=[];
%         i_temp = i_temp + nx_fine-1;
%     end
% end


% source force F
% FF = zeros((m*nt_interval+1)*(nx_fine-1),1);
% FF(1:nx_fine-1) = ux0;
% % i_temp = 1;
% % t0 = tF(1); i = i+1
% for ii = 2:m*nt_interval+1
%     FF((ii-1)*(nx_fine-1) + 1 :ii*(nx_fine-1) ,:) = A_Fo^-1*dt*f(xF(2:nx_fine),tF(ii))';
% end

% uF_seq = A_tilde\( R_rhs*F);
% Right_hand_side = R_rhs*F;
% u0 = Right_hand_side(1:nx_fine-1,1);
uF_seq =  FineSolver(A_F,m,nt_interval,nx_fine,FF);

uF_seq_mat=[];
for i = 1:length(uF_seq)/(nx_coarse-1)
    uF_seq_mat(:,i) = [uF_seq((i-1)*(nx_coarse-1)+1:(i)*(nx_coarse-1))];
end


%     Uk_0 = uF_true;
% end
% figure
% hold on
% plot(xF,[0;uF_seq(end-nx_fine+2:end);0])
% 
% plot(xF,U_exact)
% legend('U_{FineSequential}','U_{ex}')
% title('Sequentially Fine Solution and the Exact Solution')

% uF_seq = U_exact;
% return
%% ------------------------------------------------------------------
% initial coarse integration

% initial solution u0
ux0C = ux0_function(xC(2:end-1))';

% Coarse sequential propagator

Uk_0C = zeros((nt_interval+1)*(nx_coarse-1),1);
Uk_0C(1:nx_coarse-1) = ux0C;
% rhs0 = [ux0; zeros(((nt_interval+1)*(nx_coarse-1)-nx_coarse+1),1)];
%  rhs = [ux0; (A_tilde-M_tilde)*Uk_0];
% source force F
F_C = zeros((nt_interval+1)*(nx_coarse-1),1);
F_C(1:nx_coarse-1) = ux0C;
% i_temp = 1;
for i = 2:nt_interval+1
    F_C((i-1)*(nx_coarse-1) + 1 :i*(nx_coarse-1) ,:) = A_C^-1*dT*f(xC(2:nx_coarse),tC(i))';
end

% figure
% hold on
% uC = M_tilde\( F_C);
% plot(xC,[0;uC(end-nx_coarse+2:end);0])

% plot(xF,U_exact)
% legend('U_{0}','U_{ex}')
% title('Parareal Solution and the Exact Solution')

% coarse space initial
indice_R_coarse  = 0:m:m*nt_interval;
R = zeros((nx_coarse-1)*length(indice_R_coarse),(m*nt_interval+1)*(nx_coarse-1));
for i = 1:length(indice_R_coarse)
    R((i-1)*(nx_coarse-1)+1:(i)*(nx_coarse-1),indice_R_coarse(i)*(nx_coarse-1)+1:indice_R_coarse(i)*(nx_coarse-1)+(nx_coarse-1)) = Id_C;
end

% coarse correction
R'*M_tilde*R;
% return
% first subdomain  Omega_1 = {0}
R1 = zeros((nx_coarse-1),(m*nt_interval+1)*(nx_coarse-1));
R1(1:nx_coarse-1,1:nx_coarse-1) = Id_C;
A1 = R1*A_Original*R1';



% define nodes in other subdomains
N_nodes_subdomain = m;

% Number of subdomains except the first subdomain
N_subdomain = m/N_nodes_subdomain*nt_interval;%nt_fine/N_nodes_subdomain

indice_R_subdomains = zeros(N_subdomain,N_nodes_subdomain);
for i =1:size(indice_R_subdomains,1)
    indice_R_subdomains(i,:) =  (i-1)*N_nodes_subdomain+1:(i-1)*N_nodes_subdomain+N_nodes_subdomain;
end



% A_Fo = (kF-rF)*diag(ones(nx_fine-2,1),1) + (1+2*rF-c*dt)*eye(nx_fine-1) - (kF+rF)*diag(ones(nx_fine-2,1),-1);
% Id_Fo = eye(size(A_Fo));
% A_Original =  zeros(m*(nt_interval+1)*(nx_coarse-1));



% Omega_2 = {1,2} Omega_3 = {3,4}
R2_end = zeros(N_subdomain*N_nodes_subdomain*(nx_coarse-1),(N_subdomain*N_nodes_subdomain+1)*(nx_coarse-1)) ;
for i = 1:N_subdomain*N_nodes_subdomain
    R2_end((i-1)*(nx_coarse-1)+1:i*(nx_coarse-1),(i)*(nx_coarse-1)+1:(i+1)*(nx_coarse-1)) = Id_F;
end

R_sub = [];       % restriction matrices
A_sub = [];       % sub-matrices
for i = 2: N_subdomain+1
    R_sub{i} = R2_end(N_nodes_subdomain*(nx_coarse-1)*(i-2)+1:N_nodes_subdomain*(nx_coarse-1)*(i-2)+N_nodes_subdomain*(nx_coarse-1),:);
    A_sub{i} = R_sub{i}*A_Original*R_sub{i}';
end

% % define 2-level domain decomposition multiplicative preconditioner
P = R1'*A1*R1;
P1 = R1'*A1^-1*R1;
for i = 2:N_subdomain+1
    P = P + R_sub{i}'*A_sub{i}*R_sub{i};
    P1 = P1 + R_sub{i}'*A_sub{i}^-1*R_sub{i};
end
DDD = P^-1 - P1;
norm(DDD);
% return


% P_sum_inverse = R1'*A1^-1*R1;
% for i = 2:nt_interval+1
%     P_sum_inverse = P_sum_inverse + R_sub{i}'*A_sub{i}^-1*R_sub{i};
% end



%%Multiplicative
Id = eye(size(A_Original));
M_mul = P*(R'*M_tilde*R +  Id - R'*R);

%%Additive
M_add = P + (R'*M_tilde*R - R'*R);

%%% Coarse space correction
B = R'*M_tilde*R +  Id - R'*R;

% initial coarse solution
uC = M_tilde\( F_C);
uC_mat = [];


for i = 1:length(uC)/(nx_coarse-1)
    uC_mat(:,i) = [uC((i-1)*(nx_coarse-1)+1:(i)*(nx_coarse-1))];
end
% uC_mat=uC_mat';
ttrue = linspace(0,T,nt_interval*m+1);
pick = 1:m:nt_interval*m+1;
U_fine_interpolation = spline(ttrue(pick),uC_mat,ttrue);
% U_fine_interpolation = interp1(ttrue(pick),uC_mat',ttrue,'spline');

% return
U_fine = [];
for i = 1:size(U_fine_interpolation,2)
    U_fine((i-1)*(nx_coarse-1)+1:(i)*(nx_coarse-1),:) = U_fine_interpolation(:,i);
end



%% 
%     r_k = FF - A_Original*U_old;
%     U_old_temp = P\r_k;
%     DU = (R'*M_tilde*R + Id - R'*R)\U_old_temp;
 

 % residual
%     r_k = F - A_Original*U_fine;
     r_k = FF - A_Original*U_old;
    
  % first level ~ fine propagator ~ additive Schwarz in parallel
%     U_fine_temp = P\r_k;   
    for i = 1:N_subdomain+1
%         i = i + 1
        if i == 1
            x_temp_1 = R1*r_k;
            x_temp_2 = A1\x_temp_1;
            x{i} = R1'*x_temp_2;
            U_old_temp = x{i};         
        else
            x_temp_1 = R_sub{i}*r_k;
            x_temp_2 = A_sub{i}\x_temp_1;
            x{i} = R_sub{i}'*x_temp_2;  
            U_old_temp = U_old_temp + x{i};
        end            
    end

  % second level ~ coarse grid correction sequentially
    U_old_temp_1 = R*U_old_temp;
    U_old_temp_2 = M_tilde\U_old_temp_1;
    U_old_temp_3 = R'*U_old_temp_2;
    U_old_temp_4 = U_old_temp_3 + U_old_temp - R'*R*U_old_temp;
    DU =   U_old_temp_4;













