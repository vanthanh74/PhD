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
function [x,r,V,H,k,resvec,u,c] = priorgmres(A,x,r0,m,M1,M2,tol)

if(isempty(M1))
   existM1 = 0;
else
   existM1 = 1;
end
if(isempty(M2))
   existM2 = 0;
else
   existM2 = 1;
end

% Initialize V
V(:,1) = r0 / norm(r0);

for k = 1:m
   % Find w using preconditioning if available.
   if(existM2)
      w = M2 \ V(:,k);
   else
      w = V(:,k);
   end
   w = A*w;
   if(existM1)
      w = M1 \ w;
   end

   % Create next column of V and H
   for j = 1:k
      H(j,k) = V(:,j)' * w;
      w = w - H(j,k) * V(:,j);
      hh=V(:,j)' * w;
      w=w-hh*V(:,j);
   end

   H(k+1,k) = norm(w);
   V(:,k+1) = w / H(k+1,k);

end
 rhs = zeros(k+1,1);
   rhs(1) = norm(r0);
  y = H \ rhs;
  res = rhs - H * y;
  resvec=norm(res);
% Calculate solution and residual.
z=V(:,1:m) * y;
x = x + z;
r = V * res;
delta=norm(r-r0);
c=r-r0/delta;
u=delta*z;