
% Backward Euler

function [U1] = RungeKuttaFine_convdiff(delta_x,delta_t,f,u0,t0,Nx,Nt)
% clear all      ( rhs , t_current , dt_fine , U_current , N_steps )
% close all
% clc



% u_ex=@(x,t) x.*(1-x).^2*exp(-2*t);
% f=@(x,t) (4 - 8*x + 4*x.^2 - 2*x.^3).*exp(-2*t);

% Nx=20;
% Nt=5;

%  delta_x = dx;
% delta_t = dt;
% u0 = U_interpolation(n,:);
% Nx = nx_fine;
% Nt = m;
% t0 = tC(n);

global xF a b c
U0 = u0(2:end-1)';
F=zeros(Nx-1,1);


r = a/(delta_x^2);
k = b/(2*delta_x);

A = (r-k)*diag(ones(Nx-2,1),1) - (2*r-c)*eye(Nx-1) + (r+k)*diag(ones(Nx-2,1),-1);

% A = (2*diag(ones(1,Nx-1))-diag(ones(1,Nx-2),1)-diag(ones(1,Nx-2),-1))/(hx^2);
sol = [U0'];

for n=1:Nt
    
    for i=2:Nx
        F(i-1)=f(xF(i),t0);
    end
    
      k1=A*U0+F;
      k2=A*(U0+k1*delta_t/2)+F;
      k3=A*(U0+k2*delta_t/2)+F;
      k4=A*(U0+k3*delta_t)+F;
      U1=U0+(delta_t/6)*(k1+2*k2+2*k3+k4);
  
%     U1 = A\(U0 + F);   
    U0 = U1;
    sol =[sol; U1'];
    t0 = t0 + delta_t;
end

if Nt>1
    sol = [zeros(size(sol,1),1) sol zeros(size(sol,1),1)];
end
U1=sol;






