function [U1,nmv_total] = RungeKuttaFine_2DConvDiff(kappa,dx_fine,dy_fine,dt,f,u0,t0,nx_fine,ny_fine,nt,name)



global xF yF p1 p2 p3

% u0 = U_Fine0
% t0 = 0
% nt = m
% name = 'a'

U0 = zeros((nx_fine-1)*(ny_fine-1),1);
for i = 1:nx_fine-1
    for j=1:ny_fine-1
        U0((i-1)*(nx_fine-1)+j) = u0(i+1,j+1);
    end
end




r1 = kappa/(dx_fine^2);
k1 = p1/dx_fine;

r2 = kappa/(dy_fine^2);
k2 = p2/dy_fine;

% A = -r1*diag(ones((ny_fine-1)*(ny_fine-1) - 1,1),-1) ...
%     + (1+2*r1+2*r2)*eye((ny_fine-1)*(ny_fine-1)) ...
%     - r1*diag(ones((ny_fine-1)*(ny_fine-1)-1,1),1) ...
%     - r2*diag(ones((ny_fine-1)*(ny_fine-1)-(ny_fine-1),1),-(ny_fine-1)) ...
%     - r2*diag(ones((ny_fine-1)*(ny_fine-1)-(ny_fine-1),1),(ny_fine-1)) ;


A =  (r1-k1)*diag(ones((ny_fine-1)*(ny_fine-1)-1,1),1) ... 
    + (r1+k1)*diag(ones((ny_fine-1)*(ny_fine-1) - 1,1),-1) ...
    - (2*r1+2*r2 - p3)*eye((ny_fine-1)*(ny_fine-1)) ...    
    + (r2-k2)*diag(ones((ny_fine-1)*(ny_fine-1)-(ny_fine-1),1),(ny_fine-1))...
    + (r2+k2)*diag(ones((ny_fine-1)*(ny_fine-1)-(ny_fine-1),1),-(ny_fine-1))  ;


B = zeros((nx_fine-1)*(ny_fine-1),(nx_fine-1)*(ny_fine-1));
for i = 1:nx_fine-2
     B(i*(nx_fine-1),i*(nx_fine-1)+1) = -(r2-k2);
     B(i*(nx_fine-1)+1,i*(nx_fine-1)) = -(r2+k2);
end

A = A+B;
A = sparse(A);
% I = r2*eye(ny_fine-1,ny_fine-1);

% A = zeros((nx_fine-1)*(ny_fine-1),(nx_fine-1)*(ny_fine-1))

% Initialize matvec count
nmv_total = [];

sol =[u0];
% x0 = zeros(length(F),1);

F  = zeros((nx_fine-1)*(ny_fine-1),1);
% F = [];

tol = 1e-12;
% U0 = zeros(length(F),1);
for  n = 1 : nt
    %      n = n+1
   
    
    
    for i=1:nx_fine-1
        for j=1:ny_fine-1
            F((i-1)*(nx_fine-1) +j)= f(i*dx_fine,j*dy_fine,t0);
        end
    end
    
    
    %-----------------------------------------------------------------------------
    % using Krylov subspace recycling
    [L,U]=ilu(A);
    
    % Call GCRO-DR
    m = 30;%the maximum subspace dimension used by GCRODR.
    k = 10;%the number of approximate eigenvectors kept from one cycle to the next.
    
%                     [U1,resvec,~,nmv,~] = gcrodr(A,U0+F,m,k,U0,tol,[],[],name);  % nonpreconditioned
%                 [U1,resvec,~,nmv,~] = gcrodr(A,U0+F,m,k,U0,tol,L,U,name);    % preconditioned

    % Record number of matrix-vector products needed to solve this linear system
%             nmv_total = [nmv_total; nmv];
    
    % Scale convergence history to show relative residual norms
    %     resvec = resvec / norm(U0+F);
    %-----------------------------------------------------------------------------
    
    
    
    
    %-----------------------------------------------------------------------------
    % using CG
    %                 U1 = pcg(A,U0+F,tol,30,[],[],[]);             % nonpreconditioned
    %                 U1 = pcg(A,U0+F,tol,20,L,U,U0);  % preconditioned
    %-----------------------------------------------------------------------------
    
    
    % using Deflated-CG
    %         U1 = DeflCG(A,U0+F,U0);
    
    
    
    %-----------------------------------------------------------------------------
    % using GMRES
%     [U1,~,~,~,resvec] = gmres(A,U0+F,50,tol,100,[],[],U0);            % nonpreconditioned ,intial guess U0?
%     [U1,~,~,~,resvec] = gmres(A,U0+F,50,tol,100,L,U,U0);  % preconditioned

%     length(resvec)-1
%     Record number of matrix-vector products needed to solve this linear system
%     nmv_total = [nmv_total;length(resvec)-1];
    %-----------------------------------------------------------------------------
%             U1 = A\(U0 + F);

      k1=A*U0+F;
      k2=A*(U0+k1*dt/2)+F;
      k3=A*(U0+k2*dt/2)+F;
      k4=A*(U0+k3*dt)+F;
      U1=U0+(dt/6)*(k1+2*k2+2*k3+k4);
  
      
     U0 = U1;
     nmv_total = [nmv_total; 1];
     t0 = t0 + dt;
end
U_matrix = zeros((nx_fine+1),(ny_fine+1));
for i = 2:nx_fine
    for j =2:ny_fine
        U_matrix(i,j) =  U1((j-2)*(nx_fine-1)+i-1) ;
    end
end
sol =[sol; U_matrix];
U1=sol;
% nmv = length(resvec)-1;
% nmv_total = [nmv_total zeros(length(nmv_total),1)];
end