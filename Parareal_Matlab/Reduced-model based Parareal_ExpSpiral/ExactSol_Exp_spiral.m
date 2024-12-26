%   function [t,u]=Runge_Kutta4_Exp_spiral(f,tspan,u0,N)
function [t,u]=ExactSol_Exp_spiral(f,tspan,u0,N)
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
%   
%       k1=A*U0+F;
%       k2=A*(U0+k1*dt/2)+F;
%       k3=A*(U0+k2*dt/2)+F;
%       k4=A*(U0+k3*dt)+F;
%       U1=U0+(dt/6)*(k1+2*k2+2*k3+k4);
      u(n+1,1) = exp(t(n+1))*((1)*cos(100*t(n+1)) - 1*sin(100*t(n+1)));
      u(n+1,2) = exp(t(n+1))*((1)*sin(100*t(n+1)) + 1*cos(100*t(n+1)));

end


     