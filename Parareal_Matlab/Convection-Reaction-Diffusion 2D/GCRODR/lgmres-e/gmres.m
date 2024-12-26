% GMRES
%
% Solves A z = r for z, then returns x + z
% Generates Arnoldi relation A V(:,1:m) = V(:,1:m+1) H
%
% INPUT:  A      N-by-N matrix
%         X      current solution vector
%         R      N-by-1 preconditioned residual vector
%         M      number of GMRES iterations to perform
%         M1     left preconditioner for A
%         M2     right preconditioner for A
%         tol    specifies the tolerance of the method
% OUTPUT: x      updated solution vector
%         R      preconditioned residual vector
%         V      N-by-M+1 matrix containing orthogonal basis for Krylov subspace 
%         H      M+1-by-M upper Hessenburg reduction of matrix operator
%         K      number of GMRES iterations actually performed
%         RESVEC vector containing norm of residual at each iteration of GMRES
function [x,r,C,U,resvec] =gmres(A,x,r,m,tol);


% Initialize V
V(:,1) = r / norm(r);

for k = 1:m
   % Find w using preconditioning if available.
  
   w = V(:,k);
  
   w = A*w;
 

   % Create next column of V and H
   for j = 1:k
      H(j,k) = V(:,j)' * w;
      w = w - H(j,k) * V(:,j);
      hh=V(:,j)' * w;
      w=w-hh*V(:,j);
   end

   H(k+1,k) = norm(w);
   V(:,k+1) = w / H(k+1,k);
   
  % Initialize right hand side of least-squares system
   rhs = zeros(k+1,1);
   rhs(1) = norm(r);

   % Solve least squares system; Calculate residual norm
   y = H \ rhs;
   res = rhs - H * y;
   resvec(k) = norm(res);
   if resvec(k) < tol
      % Calculate solution and residual and return
      x = x + V(:,1:k) * y;
      r = V * res;
      return
   end
end
Z_p=V(:,1:m) * y;
x = x + Z_p;
r1 = V * res;
Q_p=r1-r;
[Q R]=qr(Q_p,0);
 C=Q;
U=Z_p/R;
r=r1;
















