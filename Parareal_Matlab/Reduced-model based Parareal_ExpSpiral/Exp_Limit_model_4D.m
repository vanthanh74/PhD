  function [t,u]=Exp_Limit_model_4D(tspan,u0,N)
global c
dt=(tspan(2)-tspan(1))/N;
t=(tspan(1):dt:tspan(2))';  % colon to make column vector
% u(1,:)=u0(:);
% for n=1:N
% %   u(n+1,:)=u(n,:)+dt*f(t(n),u(n,:));
%    
%       k1= dt*f(t(n),u(n,:));
%       k2= dt*f(t(n)+dt/2,u(n,:) + k1/2);
%       k3= dt*f(t(n)+dt/2,u(n,:) + k2/2);
%       k4= dt*f(t(n)+dt,u(n,:) + k3);
%       u(n+1,:) = u(n,:)+(1/6)*(k1+2*k2+2*k3+k4);
% 
% 
% end

 u(1)  = u0(1)+(1/c)*(dt*(u0(1)+2*u0(2))+sin(c*dt)*u0(3)+(1-cos(c*dt))*u0(4));
 u(2) = u0(2)+(1/c)*(-dt*(2*u0(1)+u0(2))+(cos(c*dt)-1)*u0(3)+sin(c*dt)*u0(4));
 u(3) = cos(c*dt)*u0(3)+sin(c*dt)*u0(4)+(1/c)*(-cos(c*dt)*(2*dt*u0(4))+sin(c*dt)*(2*dt*u0(3))+sin(c*dt)*(2*u0(1)+u0(2))+(1-cos(c*dt))*(u0(1)+2*u0(2)) ) ;
 u(4) = -sin(c*dt)*u0(3)+cos(c*dt)*u0(4)+(1/c)*(sin(c*dt)*(2*dt*u0(4))+cos(c*dt)*(2*dt*u0(3))+(cos(c*dt)-1)*(2*u0(1)+u0(2))+sin(c*dt)*(u0(1)+2*u0(2)) );  % Limidt model


     