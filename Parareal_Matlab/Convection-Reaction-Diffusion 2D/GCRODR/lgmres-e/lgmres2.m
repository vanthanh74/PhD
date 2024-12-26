% This code implement a variant of Lgmres method. 
% an explicit orthogonalization process is performed,
% (I-C_pC_p^{T})A is used replace A , where 
% AZ_p=QR,
% C_p=Q,
% U_p=Z_p*inv(R);
% The final relationship is A[U_p V_{m-p}]=[C_p V_{m-p+1}]*[I_p,B_p;0,H_{m-p+1}]
%                          
% so the approximate solution is by solving min||||??
% A=mmread('nos3.mtx');n=size(A,1);b=rand(n,1);x=rand(n,1);r=b-A*x;m=40;p=1;

function [iter] = lgmres2(A,x,r,m,p,tol);
tic
iter=1;
[x,r,C,U,resvec] =gmres(A,x,r,m,tol);
while resvec(iter)>tol
iter=iter+1;

% Initialize V
V(:,1) = r / norm(r);
for k = 1:m-p
   % Find w using preconditioning if available.
   w = V(:,k);
   w = A*w;

   % Apply (I-C*C') operator to Aw
   B(:,k) = C' * w;
   w = w - C * B(:,k);

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
  for i = 1:p
      d(i) = norm(U(:,i));
      U(:,i) = U(:,i) / d(i);
   end
   D = diag(1 ./ d);

   % Form large H
   H2 = zeros(m+1,m);
   H2(1:p,1:p) = D;
   H2(1:p,p+1:m) = B;
   H2(p+1:m+1,p+1:m) = H;

   % Calculate solution update
   rhs = [C V]' * r;
   y = H2 \ rhs;
   
 Z_p= [U V(:,1:m-p)] * y;
x = x + Z_p;
Q_p=[C V] * (H2 * y);
[Q R]=qr(Q_p,0);
C=Q;
U=Z_p/R;
r = r - Q_p;
resvec(iter) = norm(r);
end
toc
 semilogy(resvec);
 
