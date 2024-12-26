  function [t,u]=Runge_Kutta4_Exp_spiral(f,tspan,u0,N)

% FORWARDEULER solves system of ODEs using the Forward Euler method
%   [t,u]=ForwardEuler(f,tspan,u0,N) solves du/dt=f(t,u) with initial
%   value u0 on the time interval tspan doing N steps of Forward
%   Euler. Returns the solution in time and space in the matrix u, and
%   also the corresponding time points in the column vector t.
dt=(tspan(2)-tspan(1))/N;
t=(tspan(1):dt:tspan(2))';  % colon to make column vector
u(1,:)=u0(:);
for n=1:N
%   u(n+1,:)=u(n,:)+dt*f(t(n),u(n,:));
   
      k1= dt*f(t(n),u(n,:));
      k2= dt*f(t(n)+dt/2,u(n,:) + k1/2);
      k3= dt*f(t(n)+dt/2,u(n,:) + k2/2);
      k4= dt*f(t(n)+dt,u(n,:) + k3);
      u(n+1,:) = u(n,:)+(1/6)*(k1+2*k2+2*k3+k4);


end


     