close all
clear all
clc


% format long
global c C
c = 100;   %eps = 1/c    
C = 2;% Exp_spiral rhs
% K=25; u_init=[1 1 1 1]; T=20*pi/c; N=10;
N=16;% parareal parameters
K=16; u_init=[1 1 1 1 1 1]; T=4*2*pi/c;  
% K=15; u_init=[1 1 1 1 1/c 0]; T=4; N=4;  
DT = T/(N)
MG=1; MF= 4*80/N;

% 40*round(DT/(2*pi/c)) %min(floor(40*T/(2*pi/c)),floor(40*T/(2*pi/c)+1));    % F and G number of time steps

% figure
% plot(1:MF,ones(MF,1),'.b',1:DT/(2*pi/c):MF,1,'dr')
% dT = DT / MG;

dt = DT / MF
m = DT/dt
M = DT/MG
M/dt
% f=@(t,x) [-c*x(2)+x(1) c*x(1)+x(2)];

f=@(t,u) [u(4) ...
          u(5) ...
          u(6) ...
          -C*u(1) ...
          c*u(6)+ C/2*u(2) ...
          -c*u(5)+C/2*u(3)];
   
fG1 = @(t,u) [   u(4) ...
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
            
F=@(t0,t1,u0) SRunge_Kutta4_Exp_spiral(f,t0,t1,u0,MF); % fine solver F
% G=@(t0,t1,u0) SGRunge_Kutta4_Exp_spiral_6D_Const_PenningTrap(fG1,t0,t1,u0,MG); %  Limit model coarse solver G coarse solver G
G=@(t0,t1,u0) SGRunge_Kutta4_Exp_spiral(f,t0,t1,u0,MG); %   coarse solver G

% G=@(t0,t1,u0) SRunge_Kutta4_Filter(fG,t0,t1,u0,MG); % coarse solver G with Filter

% U0 = G(0,DT,u_init)
U=Parareal(F,G,T,u_init,N,K);                   % solve with parareal
[t,u]=Runge_Kutta4_Exp_spiral(f,[0 T],u_init,MF*N);        % fine solution
TT=0:T/N:T;                                 % coarse time mesh
% U_InitCoarse = u_init;
% for i =1:N
%     U_InitCoarse(i+1,:) = G(0,DT,U_InitCoarse(i,:));
% end

% L2error = [norm(u(1:MF:end,:)-U_InitCoarse,2)];
% LInferror = [norm(u(1:MF:end,:)-U_InitCoarse,inf)];

L2error = [norm(u(1:MF:end,:)-U{1},2)];
LInferror = [norm(u(1:MF:end,:)-U{1},inf)];
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
  L2error(k+1,:) = [norm(u(1:MF:end,:)-U{k+1},2)];
  LInferror(k+1,:) = [norm(u(1:MF:end,:)-U{k+1},inf)];
%   pause
end
legend('Fine solution','Parareal solution')
title(['Solution X  with N = ', num2str(N)])
view([90 60 30])
figure
%   plot3(u(:,4),u(:,5),u(:,6),'-b'...        % solution and         
%     ,U{k}(:,4),U{k}(:,5),U{k}(:,6),'or');
 plot(u(:,4),u(:,5),'-b'...        % solution and         
    ,U{k+1}(:,4),U{k+1}(:,5),'or');
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

% figure
% semilogy(0:K,L2error,'--^b',0:K,LInferror,'--^r','linewidth',1.5)
% legend({'L^{2}-norm','L^{\infty}-norm'},'location','northeast','fontsize',15) %northeast southwest
% title(['N = ', num2str(N),' and T = ',num2str(T)],'fontsize',15)
% xlabel('Iteration')
% ylabel('Error')

figure
semilogy(0:K,LInferror,'--^r','linewidth',1.5)
legend({'L^{\infty}-norm'},'location','northeast','fontsize',15) %northeast southwest
title(['N = ', num2str(N),' and T = ',num2str(T)],'fontsize',15)
xlabel('Iteration')
ylabel('Error')


% figure
% semilogy(0:K,L2error,'b',0:K,LInferror,'r',[2  5.3],[ 1e-6 1e-13],'k')
% legend('L^{2}error','L^{\infty}error','quadratic slope')
% title(['Convergence rate   with N = ', num2str(N)])

% figure 
% plot(u(:,4:6))
L2error
LInferror
% ylim([1e-10 1e0])
% ylim([1e-15 1e1])