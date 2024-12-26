%************************************************************************************** = 
% AUTHOR      : Michael L. Parks, University of Illinois at Urbana-Champaign
% DATE        : 3/30/04
% DESCRIPTION : Solve a sequence of 10 linear systems with GCRODR, 
%               using Krylov subspace recycling (on systems 2-10)
% NOTES       : Developed by Michael L. Parks and Eric de Sturler, University of Illinois
%               Contact: parks@cse.uiuc.edu,  http://www.cse.uiuc.edu/~parks 
%                        sturler@cs.uiuc.edu, http://www-faculty.cs.uiuc.edu/~sturler/
%               See UIUC CS tech report UIUCDCS-R-2004-2421
% ACK         : Example linear systems arise from finite element analysis of crack propogation.
%               Problem courtesy of Philippe H. Geubelle (AAE-UIUC, geubelle@uiuc.edu) 
%               and Spandan Maiti (AAE-UIUC)
%**************************************************************************************

% Be sure to "clear all" so that no subspace is recycled from previous calls to driver 

close all;clc;clear all;
load('dataAll3.mat');
% Put GCRODR in path
% addpath('gcrodr');

% Initialize count
count = 0;

% Initialize matvec count
nmv_total = [];

A2=Kall{2};
 
%[L U]=luinc(A,'0');
% Solve 10 sequential linear systems
tic  

x0 = zeros(length(F),1);
for i = 1:10
    
  % disp(sprintf('\nNow solving linear system #%i\n',i))
 % b=I(:,i);
   % Load system from disk
   
   Ai=Kall{i};
    
    A=Ai(1:end-6,1:end-6);
  [L U]=ilu(Ai);
 %   D=[A-tril(A,-2)-triu(A,2),b;zeros(size(b')),c-b'*diag(diag(A))*b];
     
   %A=D\A;
  % x0 = zeros(length(b),1);
   % Call GCRO-DR
   [x0,resvec,r,nmv,relres] = gcrodr(Ai,F,20,10,x0,1e-10,L,U);

   % Record number of matrix-vector products needed to solve this linear system
   nmv_total = [nmv_total, nmv];
   
   % Scale convergence history to show relative residual norms
   resvec = resvec / norm(F);

   % Plot convergence history
   u = count : count+length(resvec)-1;
   semilogy(u,resvec,'-.b');
   axis([0 4000 1e-10 1]);
   title('Convergence for 10 consecutive linear systems')   
   xlabel('Total number of matrix-vector products')
   ylabel('Relative residual norm')
   hold on;
toc
   count = count + length(resvec);
   
end

% Output number of matrix-vector products needed to solve each linear system
for i = 1:10
   disp(sprintf('Linear system %2i required %i matrix-vector products.',i,nmv_total(i)))
end