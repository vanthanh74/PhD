close all
clear all
clc

% IVP: u_t = kappa*(u_xx + u_yy) - 2*p1*u_x - 2*p2*u_y + p3*u + f    in [0,L]x[0,L] x (0,T),
%      u(0,y,t) = u(L,y,t)  = 0                 0 <= y <= L, t >= 0
%      u(x,0,t) = u(x,L,t)  = 0                 0 <= x <= L, t >= 0
%      u(x,y,0) = u0(x,y)                   (x,y) in [0,L]x[0,L]

% INPUT PARAMETERS:
%   * nt_coarse: number of coarse intervals in time
%   * nx_coarse, ny_coarse: number of coarse intervals in space
%   * nt_fine: number of fine intervals in time
%   * nx_fine, ny_fine: number of fine intervals in space
%   * m: number of fine time steps on each coarse time step 
%   * K: number of Parareal iterations (2 <= K <= n_coarse+1)
%        (--> K sequential coarse sweeps and K-1 parallel fine sweeps)
% OUTPUT:
%   * LInfinityErrorNorm

% Put GCRODR in path
addpath('gcrodr');

global T L kappa p1 p2 p3 

 T = 1; %  Intervall (0,T)
 L = 1; %  Omega=(0,L)
 
 kappa = 1;
 p1 = 0.05;
 p2 = 0.05;
 p3 = 1;

nx_coarse = 10;
ny_coarse = 10; 
nt_coarse = 20;

nx_fine = 20;
ny_fine = 20;
nt_fine = 200;


dx_fine = L/nx_fine;   % fine spatial discretization steps
dy_fine = L/ny_fine;   % fine spatial discretization steps
dt = T/nt_fine;   % fine temporal discretization steps % dx_fine^2/4 for Runge-Kutta


dX_coarse = L/nx_coarse; % coarse spatial discretization steps
dY_coarse = L/ny_coarse; % coarse spatial discretization steps
dT = T/nt_coarse; % coarse temporal discretization steps

m = round(dT/dt)  % number of fine time steps on each coarse time step 


K= 20   % number of Parareal iterations

% LInfNormError = parareal_1d_heat_testAlpha(L,T,dT,dX,dt,dx,nt_coarse,nx_coarse,nt_fine,nx_fine,m,K)
[nmv_total_matrix_FineSolver,nmv_total_matrix_FineSolverInParareal,TimeConsumingFineSolver,TimeConsumingParareal,total_nmv_perIterParareal,total_nmv] = parareal_2d_convection_diffusion(L,T,kappa,dT,dX_coarse,dY_coarse,...
                                 dt,dx_fine,dy_fine,nt_coarse,nx_coarse,ny_coarse,...
                                 nt_fine,nx_fine,ny_fine,m,K)
