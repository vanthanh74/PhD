close all
clear all
clc

%  X1'(t) = V1(t)
%  X2'(t) = V2(t)
%  X3'(t) = V3(t)
%  V1'(t) = 1/eps*(-X1*V3/sqrt(X1^2+X2^2)) + V2      ,   0<eps<<1, eps=0.01
%  V2'(t) = 1/eps*(-X2*V3/sqrt(X1^2+X2^2)) - V1
%  V3'(t) = 1/eps*((X1*V1 + X2*V2)/sqrt(X1^2+X2^2))
%  X = (X1,X2,X3)'  ,  V = (V1,V2,V3)'



% format long
global c C
c = 100;   %eps = 1/c   
C = 2;% Exp_spiral rhs
% K=25; u_init=[1 1 1 1]; T=20*pi/c; N=10;
N=2;% parareal parameters
iter_refine = 7;
% K=20; u_init=[1 1 1 1 1 1 0 0 0 0 0 0 ]; T=10;
 K=20; u_init=[1 1 1 1 1 1  ]; T=5;
% K=15; u_init=[1 1 1 1 1/c 0]; T=4; N=4;
DT = T/(N)
MG=1;
% 40*round(DT/(2*pi/c)) %min(floor(40*T/(2*pi/c)),floor(40*T/(2*pi/c)+1));    % F and G number of time steps

% figure
% plot(1:MF,ones(MF,1),'.b',1:DT/(2*pi/c):MF,1,'dr')
% dT = DT / MG;

% f=@(t,x) [-c*x(2)+x(1) c*x(1)+x(2)];

f=@(t,u) [u(4) ...
          u(5) ...
          u(6) ...
          -C*u(1) ...
          c*u(6)+ C/2*u(2) ...
          -c*u(5)+C/2*u(3)];
   
fG = @(t,u) [   u(4) ...
                0 ...
                0 ...
                -C*u(1) ...
               0 ...
                0];




%
% fG2 = @(t,u,u0) [   u(4) ...
%                -u0(3)+u0(1) ...
%                 u0(2)-u0(1) ...
%                 -u(1) ...
%                u0(6)+u0(1) ...
%                 -u0(5)-u0(4)];
L2error = zeros(K,iter_refine);
LInferror = zeros(K,iter_refine);
N_vec = [N];
for ii = 1:iter_refine
    MF= 12800/N;
    dt = DT / MF
    m = DT/dt
    
    F=@(t0,t1,u0) SRunge_Kutta4_Exp_spiral(f,t0,t1,u0,MF); % fine solver F
    G=@(t0,t1,u0) SGRunge_Kutta4_Exp_spiral_6D_Const_PenningTrap(fG,t0,t1,u0,MG); %  Limit model coarse solver G coarse solver G
%     G=@(t0,t1,u0) SGRunge_Kutta4_Exp_spiral(f,t0,t1,u0,MG); %   coarse solver G
    
    % G=@(t0,t1,u0) SRunge_Kutta4_Filter(fG,t0,t1,u0,MG); % coarse solver G with Filter
    
    % U0 = G(0,DT,u_init)
    U=Parareal(F,G,T,u_init,N,K);                   % solve with parareal
    [t,u]=Runge_Kutta4_Exp_spiral(f,[0 T],u_init,MF*N);        % fine solution
    TT=0:T/N:T;                                 % coarse time mesh
    U_InitCoarse = u_init;
    for i =1:N
        U_InitCoarse(i+1,:) = G(0,DT,U_InitCoarse(i,:));
    end
    
    % L2error = [norm(u(1:MF:end,:)-U_InitCoarse,2)];
    % LInferror = [norm(u(1:MF:end,:)-U_InitCoarse,inf)];
    
    %     L2error = [norm(u(1:MF:end,:)-U{1},2)];
    %     LInferror = [norm(u(1:MF:end,:)-U{1},inf)];
    L2error(1,ii) = [norm(u(1:MF:end,:)-U{1},2)];
    LInferror(1,ii) = [norm(u(1:MF:end,:)-U{1},inf)];
    figure
    for k=1:K                                 % plot fine
        %     subplot(1,2,1)
        plot3(u(:,1),u(:,2),u(:,3),'-b'...        % solution and
            ,U{k}(:,1),U{k}(:,2),U{k}(:,3),'or');    % parareal iterate
        
        xlabel('x'); ylabel('y');  zlabel('z');
        %   axis([-25 25 -25 25 ]);
        %   subplot(1,2,2)
        %   plot3(u(:,4),u(:,5),u(:,6),'-b'...        % solution and
        %     ,U{k}(:,4),U{k}(:,5),U{k}(:,6),'or');
        %
        %   xlabel('x'); ylabel('y');  zlabel('z');
        grid on
        L2error(k+1,ii) = norm(u(1:MF:end,:)-U{k+1},2);
        LInferror(k+1,ii) = norm(u(1:MF:end,:)-U{k+1},inf);
        
        %   pause
    end
    legend('Fine solution','Parareal solution')
    title(['Solution X  with N = ', num2str(N)])
    view([90 60 30])
    figure
    plot3(u(:,4),u(:,5),u(:,6),'-b'...        % solution and
        ,U{k}(:,4),U{k}(:,5),U{k}(:,6),'or');
    legend('Fine solution','Parareal solution')
    title(['Solution V  with N = ', num2str(N)])
    % figure
    % plot(t,u(:,1),'r',t(1:MF:end),U{end}(:,1),'ob')
    % legend('Fine solution','Parareal solution')
    %
    % figure
    % plot(t,u(:,2),'r',t(1:MF:end),U{end}(:,2),'ob')
    % legend('Fine solution','Parareal solution')
    %
    % figure
    % plot(t,u(:,3),'r',t(1:MF:end),U{end}(:,1),'ob')
    % legend('Fine solution','Parareal solution')
    N=N*2;
    N_vec = [N_vec;N];
end
N_vec(end) =[];

% figure
% semilogy(0:K,L2error,'b',0:K,LInferror,'r')
% legend('L^{2}error','L^{\infty}error')
% title(['Convergence rate   with N = ', num2str(N)])

leg = cell(1,iter_refine);
figure
hold on
for i = 1:iter_refine
    semilogy(0:K,LInferror(:,i),'-o')  ;
    leg{i}=['L^{\infty}error with N = ',num2str(N_vec(i))];
end
set(gca,'yscale','log')
legend(leg)
title(['Convergence rate  with T = ', num2str(T)])

% figure
% semilogy(0:K,L2error,'b',0:K,LInferror,'r',[2  5.3],[ 1e-6 1e-13],'k')
% legend('L^{2}error','L^{\infty}error','quadratic slope')
% title(['Convergence rate   with N = ', num2str(N)])

L2error
LInferror
% ylim([1e-12 1e4])