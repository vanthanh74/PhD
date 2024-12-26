close all
clear all
clc

%  X'(t) = V(t)        in R^2
%  V'(t) = 1/eps*V_per(t) + E(X)     in R^2  ,   0<eps<<1, eps=0.01         
%  X = (x,y)'  ,  V = (v_x,v_y)'  ,   V_per = (v_y,-v_x)'  and E(X) = (2x+y,x+2y)'   



% format long
global c 
c  = 100;   %eps                              % Exp_spiral rhs
ep = 1/c;
T=10; 
N=128; 

a_ep = sqrt((1-4*ep^2-sqrt(1-8*ep^2+4*ep^4))/(2*ep^2));
u_ep = 2 + a_ep^2;
w_ep = 1 + a_ep^2/ep^2;


% K=25; u_init=[1 1 1 1]; T=20*pi/c; N=10;             % parareal parameters
K=N; 
u_init=[1 1 1 1];  
% u_init=[1 0 ep -2*ep*u_ep/(2-ep^2)+ep^3*u_ep/(2-ep^2)]; 
% u_init=[1 -u_ep/w_ep ep*w_ep/u_ep -2*ep*u_ep/(2-ep^2)+ep*w_ep/u_ep+ep^3*u_ep/(2-ep^2)]; 
 




DT = T/(N)
MG=1; MF= 12800/N
% 40*round(DT/(2*pi/c)) %min(floor(40*T/(2*pi/c)),floor(40*T/(2*pi/c)+1));    % F and G number of time steps

figure
plot(1:MF,ones(MF,1),'.b',1:DT/(2*pi/c):MF,1,'dr')
% dT = DT / MG;
dt = DT / MF

% f=@(t,x) [-c*x(2)+x(1) c*x(1)+x(2)];

f=@(t,u) [u(3) u(4) c*u(4)+2*u(1)+u(2) -c*u(3)+u(1)+2*u(2)];

% KCos = @(t) (1+cos(2*pi*(t-0.5)))*(t>=0 & t<=1);

% fG = @(t0,t,x)   [-c*x(2)+x(1)*KCos((t-t0)/DT) c*x(1)+KCos((t-t0)/DT)*x(2)];

% fG2 =  @(t,u) [u_init(1)+(1/c)*(t*(u_init(1)+2*u_init(2))+sin(c*t)*u_init(3)+(1-cos(c*t))*u_init(4))...
%                  u_init(2)+(1/c)*(-t*(2*u_init(1)+u_init(2))+(cos(c*t)-1)*u_init(3)+sin(c*t)*u_init(4))...
%                  cos(c*t)*u_init(3)+sin(c*t)*u_init(4)+(1/c)*(-cos(c*t)*(2*t*u_init(4))+sin(c*t)*(2*t*u_init(3))+sin(c*t)*(2*u_init(1)+u_init(2))+(1-cos(c*t))*(u_init(1)+2*u_init(2)) ) ...
%                  -sin(c*t)*u_init(3)+cos(c*t)*u_init(4)+(1/c)*(sin(c*t)*(2*t*u_init(4))+cos(c*t)*(2*t*u_init(3))+(cos(c*t)-1)*(2*u_init(1)+u_init(2))+sin(c*t)*(u_init(1)+2*u_init(2)) )];  % Limit model

F=@(t0,t1,u0) SRunge_Kutta4_Exp_spiral(f,t0,t1,u0,MF); % fine solver F
G=@(t0,t1,u0) SRunge_Kutta4_Exp_spiral(f,t0,t1,u0,MG); % coarse solver G
% G=@(t0,t1,u0) SExp_Limit_model_4D(t0,t1,u0,MG); % Limit model coarse solver G

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

L2error = [norm(u(1:MF:end,:)-U{1},2)];
LInferror = [norm(u(1:MF:end,:)-U{1},inf)];
figure
for k=1:K                                 % plot fine    
  plot(u(:,1),u(:,2),'-b'...        % solution and         
    ,U{k}(:,1),U{k}(:,2),'or');    % parareal iterate
%   axis([-5 5 -5 5 ]); 
  xlabel('x'); ylabel('y'); 
  grid on
  L2error(k+1,:) = [norm(u(1:MF:end,:)-U{k+1},2)];
  LInferror(k+1,:) = [norm(u(1:MF:end,:)-U{k+1},inf)];
%   pause
end
legend('Fine solution','Parareal solution')
title(['Solution X  with N = ', num2str(N)])

figure
plot(u(:,3),u(:,4),'r',U{end}(:,3),U{end}(:,4),'ob')
legend('Fine solution','Parareal solution')
title(['Solution V  with N = ', num2str(N)])


% figure
% plot(t,u,'r',t(1:MF:end),U{end},'ob')
% legend('Fine solution','Parareal solution')
% title(['Final Solution  with N = ', num2str(N)])


figure
semilogy(0:K,L2error,'b',0:K,LInferror,'r')
legend('L^{2}error','L^{\infty}error')
title(['Convergence rate   with N = ', num2str(N)])

% figure
% semilogy(0:K,L2error,'b',0:K,LInferror,'r',[1  3.5 ],[ 1e-6 1e-14],'k')
% legend('L^{2}error','L^{\infty}error','quadratic slope')
% title(['Convergence rate   with N = ', num2str(N)])

% figure
% plot(log(0:K),log(L2error),'b',log(0:K),2*log(L2error)+35,'r',...
%     log(0:K),1.5*log(L2error)+17,'g',log(0:K),1.3*log(L2error)+10,'m',...
%     log(0:K),1*log(L2error)+8,'k')
% legend('L^{2}error','2log(L2error)','1.5log(L2error)','1.3log(L2error)','log(L2error)')

% figure
% plot(log(0:K),-log(L2error),'b',log(0:K),-log(LInferror),'r',log(0:K),log(0:K).^4+5,'k')
% legend('L^{2}error','L^{\infty}error','log(L2error)^4')
% xlabel('log(0:K)')
% ylabel('log(Error)')

L2error

LInferror

Error.L2error = L2error;
Error.LInferror = LInferror;

save('Error.mat','-struct','Error')
% ylim([1e-12 1e4])