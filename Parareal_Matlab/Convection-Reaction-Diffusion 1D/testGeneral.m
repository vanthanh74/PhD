close all
clear 
clc

% IVP: u_t = a*u_xx - b*u_x + c*u + f               in (0,L) x (0,T),
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

%--------------------------------------
% solution x(L-x)^2

a = 0.01;
b = 1;
c = 100;

nt_interval = 20; % number of time intervals

nx_coarse = 20; % number of coarse spatial intervals
nt_coarse = 20;% 300 %20    % number of coarse time intervals


nx_fine = 	20; % number of fine spatial intervals
nt_fine = 40;% 3000 %200 . % number of fine time intervals (nt_coarsex2 )

%--------------------------------------
% solution sin

% a = 0.5;
% b = 0.0025;
% c = 0;
% 
% nt_interval = 20; % number of time intervals
% 
% nx_coarse = 20;
% nt_coarse = 20;% 300 %20    % number of coarse time intervals
% 
% 
% nx_fine = 	20;
% nt_fine = 400;% 3000 %200 . % number of fine time intervals




% solver ='BackwardEuler';
solver ='BackwardEuler';

dx = L/nx_fine;   % fine spatial discretization steps
if strcmp(solver,'BackwardEuler')==1
    dt = T/nt_fine;   % fine temporal discretization steps . % dx^2/2 for Runge-Kutta
elseif strcmp(solver,'Runge-Kutta4')==1
    dt = dx^2/6;
end

DT = T/nt_interval;  % reference temporal discretization steps

dX = L/nx_coarse; % coarse spatial discretization steps
dT = T/nt_coarse; % coarse temporal discretization steps

M = DT/dT  % number of coarse time steps on each time interval 
m = DT/dt  % number of fine time steps on each time interval 

r = dT/dX^2;
r1 = dt/dx^2;

K = 20    % number of Parareal iterations

                                                                                                     

% LInfNormError = parareal_1d_heat_testAlpha(L,T,dT,dX,dt,dx,nt_coarse,nx_coarse,nt_fine,nx_fine,m,K)
 parareal_1d_conv_diff_General(L,T,dT,dX,dt,dx,nt_interval,nt_coarse,nx_coarse,nt_fine,nx_fine,M,m,K,solver)
