close all
clear all
clc



%(x1,x2)

% format long
global c 
c = 100;                                 % Exp_spiral rhs
K=40; u_init=[1;1]; T=80*pi/c; N=4;             % parareal parameters
DT = T/(N)
MG=1; MF= 40*round(DT/(2*pi/c)) %min(floor(40*T/(2*pi/c)),floor(40*T/(2*pi/c)+1));    % F and G number of time steps

figure
plot(1:MF,ones(MF,1),'.b',1:DT/(2*pi/c):MF,1,'dr')
% dT = DT / MG;
% dt = DT / MF;

f=@(t,x) [-c*x(2)+x(1) c*x(1)+x(2)];

KCos = @(t) (1+cos(2*pi*(t-0.5)))*(t>=0 & t<=1);

% fG = @(t0,t,x)   [-c*x(2)+x(1)*KCos((t-t0)/DT) c*x(1)+KCos((t-t0)/DT)*x(2)];
fG2 =  @(t,x) [x(1) x(2)];  % Limit model

F=@(t0,t1,u0) SRunge_Kutta4_Exp_spiral(f,t0,t1,u0,MF); % fine solver F
% G=@(t0,t1,u0) SRunge_Kutta4_Exp_spiral(f,t0,t1,u0,MG); % coarse solver G
G=@(t0,t1,u0) SGRunge_Kutta4_Exp_spiral(fG2,t0,t1,u0,MG); % Limit model coarse solver G



% G=@(t0,t1,u0) SRunge_Kutta4_Filter(fG,t0,t1,u0,MG); % coarse solver G with Filter

% U0 = G(0,DT,u_init)
U=Parareal(F,G,T,u_init,N,K);                   % solve with parareal
[t,u]=Runge_Kutta4_Exp_spiral(f,[0 T],u_init,MF*N);        % fine solution
TT=0:T/N:T;                                 % coarse time mesh
U_InitCoarse = u_init';
for i =1:N
    U_InitCoarse(i+1,:) = G(0,DT,U_InitCoarse(i,:));
end

% L2error = [norm(u(1:MF:end,:)-U_InitCoarse,2)];
% LInferror = [norm(u(1:MF:end,:)-U_InitCoarse,inf)];

L2error = [norm(u(1:MF:end,:)-U{1},2)];
LInferror = [norm(u(1:MF:end,:)-U{1},inf)];
for k=1:K                                 % plot fine    
  plot(u(:,1),u(:,2),'-b'...        % solution and         
    ,U{k}(:,1),U{k}(:,2),'or');    % parareal iterate
%   axis([-25 25 -25 25 ]); 
  xlabel('x'); ylabel('y'); 
  grid on
  L2error(k+1,:) = [norm(u(1:MF:end,:)-U{k+1},2)];
  LInferror(k+1,:) = [norm(u(1:MF:end,:)-U{k+1},inf)];
%   pause
end
legend('Fine solution','Parareal solution')
title(['Solution  with N = ', num2str(N)])

figure
plot(t,u,'r',t(1:MF:end),U{end},'ob')
legend('Fine solution','Parareal solution')
title(['Final Solution  with N = ', num2str(N)])

figure
semilogy(0:K,L2error,0:K,LInferror)
legend('L^{2}error','L^{\infty}error')
title(['Convergence rate with N = ', num2str(N)])



% ylim([1e-12 1e4])