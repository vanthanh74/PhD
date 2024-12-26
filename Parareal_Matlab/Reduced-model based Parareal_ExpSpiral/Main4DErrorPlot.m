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
global c 
c = 100;   %eps = 1/c                              % Exp_spiral rhs

K=15; u_init=[1 1 1 1]; T=10; N=20; 
Ntimes = 5;   % number of time refining N

MG=1; 


DTmat = [];
dtmax = [];
% figure
% plot(1:MF,ones(MF,1),'.b',1:DT/(2*pi/c):MF,1,'dr')
% dT = DT / MG;
% dt = DT / MF

% f=@(t,x) [-c*x(2)+x(1) c*x(1)+x(2)];

f=@(t,u) [u(3) u(4) c*u(4)+2*u(1)+u(2) -c*u(3)+u(1)+2*u(2)];
      
% fG = @(t,u) [  (u(2)^2*u(4)-u(1)*u(2)*u(5))*(1/(u(1)^2+u(2)^2)) ...
%                (-u(1)*u(2)*u(4)+u(1)^2*u(5))*(1/(u(1)^2+u(2)^2)) ...
%                 0 ...
%                 (u(5)*u(4)*u(2)-u(5)^2*u(1))*(1/(u(1)^2+u(2)^2)) ...
%                 (u(4)*u(5)*u(1)-u(4)^2*u(2))*(1/(u(1)^2+u(2)^2)) ...
%                 0];


G=@(t0,t1,u0) SExp_Limit_model_4D(t0,t1,u0,MG);
% G=@(t0,t1,u0) SExp_Limit_model_6D(t0,t1,u0,MG); % Limit model coarse solver G

% G=@(t0,t1,u0) SRunge_Kutta4_Filter(fG,t0,t1,u0,MG); % coarse solver G with Filter

cmap = hsv(Ntimes);
ErrorMax = [];
% str = {strcat('N = ', num2str(N))}
str = {};
for inum = 1:Ntimes
   
     DT = T/(N);
       MF= 100*160/N
     dt = DT / MF; 
    DT_dt=DT/dt
     DTmat(inum,1) = DT;
     dtmat(inum,1) = dt;
  
     [N MF MF*N]
    F=@(t0,t1,u0) SRunge_Kutta4_Exp_spiral(f,t0,t1,u0,MF); % fine solver F

% U0 = G(0,DT,u_init)
U=Parareal(F,G,T,u_init,N,K);                   % solve with parareal
[t,u]=Runge_Kutta4_Exp_spiral(f,[0 T],u_init,MF*N);        % fine solution
TT=0:T/N:T;                                 % coarse time mesh


% U_InitCoarse = u_init;
% for i =1:N
%     U_InitCoarse(i+1,:) = G((i-1)*DT,i*DT,U_InitCoarse(i,:));
% end

% L2error = [norm(u(1:MF:end,:)-U_InitCoarse,2)];
% LInferror = [norm(u(1:MF:end,:)-U_InitCoarse,inf)];

LInferror = [norm(u(1:MF:end,:)-U{1},inf)];

% figure
for k=1:K                                 % plot fine   
%     subplot(1,2,1)
%   plot3(u(:,1),u(:,2),u(:,3),'-b'...        % solution and         
%     ,U{k}(:,1),U{k}(:,2),U{k}(:,3),'or');    % parareal iterate
% 
% xlabel('x'); ylabel('y');  zlabel('z');
%   axis([-25 25 -25 25 ]); 
%   subplot(1,2,2)
%   plot3(u(:,4),u(:,5),u(:,6),'-b'...        % solution and         
%     ,U{k}(:,4),U{k}(:,5),U{k}(:,6),'or');
% 
%   xlabel('x'); ylabel('y');  zlabel('z');
%   grid on
%   L2error(k+1,:) = [norm(u(1:MF:end,:)-U{k+1},2)];
  LInferror(k+1,:) = [norm(u(1:MF:end,:)-U{k+1},inf)];
%   pause
end



% legend('Fine solution','Parareal solution')
% title(['Solution X  with N = ', num2str(N)])
% view([90 60 30])
% figure
%   plot3(u(:,4),u(:,5),u(:,6),'-b'...        % solution and         
%     ,U{k}(:,4),U{k}(:,5),U{k}(:,6),'or');
% legend('Fine solution','Parareal solution')
% title(['Solution V  with N = ', num2str(N)])


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

N = N*2;

ErrorMax = [ErrorMax LInferror];
% if Ntimes>3
%     break
% end

str = [str,strcat('N = ', num2str(N/2))]

end



figure
% semilogy([2  5.3],[ 1e-6 1e-13],'k')
% hold on

% Kmatrix =  (0:K).*ones(Ntimes,K+1);
for i=1:Ntimes
    semilogy(0:K,ErrorMax(:,i),'color',cmap(i,:),'Marker','o')
    hold on
end
legend(str{:}) 
   title('Convergence rate ')
   ErrorMax
%    Row={};
%    for i=1:K+1
%       Row{i} = ['K= ',num2str(i-1)] ;
%    end
% 
% T=table(ErrorMax,'RowNames',Row)

% ylim([1e-12 1e4])
DTmat
dtmat