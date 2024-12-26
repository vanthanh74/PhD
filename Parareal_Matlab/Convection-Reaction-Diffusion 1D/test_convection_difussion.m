close all
clear 
clc

% IVP: u_t = a*u_xx - b*u_x + c*u  + f               in (0,L) x (0,T),
%      u(x,0) = u0(x)                    in (0,L),
%      u(0,t) = u(L,t) = 0                 t in (0,T)  


% INPUT PARAMETERS:
%   * nt_coarse: number of coarse intervals in time
%   * nx_coarse: number of coarse intervals in space
%   * nt_fine: number of fine intervals in time
%   * nx_fine: number of fine intervals  in space
%   * m: number of fine time steps on each coarse time step 
%   * K: number of Parareal iterations (2 <= K <= n_coarse+1)
%        (--> K sequential coarse sweeps and K-1 parallel fine sweeps)
% OUTPUT:
%   * LInfinityErrorNorm

global a b c

T = 1; %  Intervall (0,T)
L = 1; % omega=(0,L)

a = 3;
b = 0.0005;
c = 1;

nx_coarse = 10;
nt_coarse = 20;% 300 %20

nx_fine = 	40;
nt_fine = 300;% 3000 %200


% solver ='BackwardEuler';
solver ='BackwardEuler';

dx = L/nx_fine;   % fine spatial discretization steps
if strcmp(solver,'BackwardEuler')==1
    dt = T/nt_fine;   % fine temporal discretization steps . % dx^2/2 for Runge-Kutta
elseif strcmp(solver,'Runge-Kutta4')==1
    dt = dx^2/6;
end

dX = L/nx_coarse; % coarse spatial discretization steps
dT = T/nt_coarse; % coarse temporal discretization steps

m = dT/dt  % number of fine time steps on each coarse time step 

r = dT/dX^2;
r1 = dt/dx^2;

K = 20    % number of Parareal iterations



% LInfNormError = parareal_1d_heat_testAlpha(L,T,dT,dX,dt,dx,nt_coarse,nx_coarse,nt_fine,nx_fine,m,K)
 parareal_1d_conv_diff(L,T,dT,dX,dt,dx,nt_coarse,nx_coarse,nt_fine,nx_fine,m,K,solver)
