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

K=15; u_init=[1 1 1 1 1/c 0]; T=4; N=10; 
Ntimes = 4 ;   % number of time refining N

MG=1;  %min(floor(40*T/(2*pi/c)),floor(40*T/(2*pi/c)+1));    % F and G number of time steps



f=@(t,u) [u(4) ...
          u(5) ...
          u(6) ...
          c*(-u(1)*u(6))/(sqrt(u(1)^2+u(2)^2))+u(5) ...
          c*(-u(2)*u(6))/(sqrt(u(1)^2+u(2)^2))-u(4) ...
          c*(u(1)*u(4) + u(2)*u(5))/(sqrt(u(1)^2+u(2)^2))];
      
fG = @(t,u) [  (u(2)^2*u(4)-u(1)*u(2)*u(5))*(1/(u(1)^2+u(2)^2)) ...
               (-u(1)*u(2)*u(4)+u(1)^2*u(5))*(1/(u(1)^2+u(2)^2)) ...
                0 ...
                (u(5)*u(4)*u(2)-u(5)^2*u(1))*(1/(u(1)^2+u(2)^2)) ...
                (u(4)*u(5)*u(1)-u(4)^2*u(2))*(1/(u(1)^2+u(2)^2)) ...
                0];
            
            
G=@(t0,t1,u0) SGRunge_Kutta4_Exp_spiral_6D(fG,t0,t1,u0,MG);
 
DT = T/N;
v_mat = [];
w_mat = [];
coef_mat =[];
TT = 0:DT:T
figure
hold on
for i = 1:length(TT)-1    
    v = rand(1,6)+1
    w = rand(1,6)+1
    v_mat = [v_mat; v]
    w_mat = [w_mat; w]
    coef = norm(G((i),TT(i+1),v)-G(TT(i),TT(i+1),w),inf)/norm(v-w,inf)  
    coef_mat = [coef_mat;coef]
    plot(TT(i),coef,'ro')
end
    