%% Heat equation del(u)/del(t) = del^2(u)/del(x^2) + f(x,t)
% Forward Euler

function [U1] = ForwardEulerFine_heat(delta_x,delta_t,f,u0,t0,Nx,Nt)
% clear all      ( rhs , t_current , dt_fine , U_current , N_steps )
% close all
% clc



% u_ex=@(x,t) x.*(1-x).^2*exp(-2*t);
% f=@(x,t) (4 - 8*x + 4*x.^2 - 2*x.^3).*exp(-2*t);

% Nx=40;
% Nt=5;



% delta_x = dx;
% delta_t = dt;
% u0 = U_interpolation(n,:);
% Nx = nx_fine;
% Nt = m;
% t0 = tC(n);


global xF
U0 = u0(2:end-1)';
F=zeros(Nx-1,1);

r = delta_t/(delta_x^2);
A = r*diag(ones(Nx-2,1),-1) + (1-2*r)*eye(Nx-1) + r*diag(ones(Nx-2,1),1);


sol = [U0'];

for n=1:Nt

    % n= n+1
%     for i=2:Nx
%         U1(i)=r*U0(i-1)+(1-2*r)*U0(i)+r*U0(i+1) + delta_t*f(x(i),t(n));
%     end

    for i=2:Nx
        F(i-1)=delta_t*f(xF(i),t0);
    end
    
    U1 = A*U0 + F;
    U0 = U1;
    t0 = t0 + delta_t;    
    sol =[sol; U1'];
    
end

if Nt>1
    sol = [zeros(size(sol,1),1) sol zeros(size(sol,1),1)];
end
U1=sol;





    
    