
% Backward Euler

function [U1] = BackwardEuler_convdiff(x,delta_x,delta_t,f,u0,t0,Nx,Nt)
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

global  a b c
U0 = u0(2:end-1)';
F=zeros(Nx-1,1);


r = a*delta_t/(delta_x^2);
k = b*delta_t/(2*delta_x);

A = (k-r)*diag(ones(Nx-2,1),1) + (1+2*r-c*delta_t)*eye(Nx-1) - (k+r)*diag(ones(Nx-2,1),-1);
sol = [U0'];

for n=1:Nt
    t0 = t0 + delta_t;
    for i=2:Nx
        F(i-1)=delta_t*f(x(i),t0);
    end
    
    U1 = A\(U0 + F);   
    U0 = U1;
    sol =[sol; U1'];
    
end

if Nt>=1
    sol = [zeros(size(sol,1),1) sol zeros(size(sol,1),1)];
end
U1=sol;






