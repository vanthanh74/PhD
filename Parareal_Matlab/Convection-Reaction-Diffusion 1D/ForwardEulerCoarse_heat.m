%% Heat equation del(u)/del(t) = del^2(u)/del(x^2) + f(x,t)
% Forward Euler

function [U1] = ForwardEulerCoarse_heat(delta_x,delta_t,f,u0,t0,Nx)
% clear all      ( rhs , t_current , dt_fine , U_current , N_steps )
% close all
% clc



% u_ex=@(x,t) x.*(1-x).^2*exp(-2*t);
% f=@(x,t) (4 - 8*x + 4*x.^2 - 2*x.^3).*exp(-2*t);

% Nx=20;
% Nt=5;

% delta_x = dX;
% delta_t = dX^2;
% u0 = ux0_restriction;
% Nx = nx_coarse;
% Nt = nt_coarse;

global xC
U0 = u0(2:end-1)';
F=zeros(Nx-1,1);


r = delta_t/(delta_x^2);
A = r*diag(ones(Nx-2,1),-1) + (1-2*r)*eye(Nx-1) + r*diag(ones(Nx-2,1),1);
U1=zeros(Nx+1,1);

for i=2:Nx
    F(i-1)=delta_t*f(xC(i),t0);
end

U1 = A*U0 + F;
U1 =[0;U1;0];







