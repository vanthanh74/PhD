% getHarmVecs2     For use with GCRODR
%
% Determines harmonic Ritz vectors using matrices computed from
% GMRES iteration. 
% 
% INPUT:  M        dimension of upper Hessenburg matrix H
%         K        select and return basis for space spanned by K harmonic 
%                  Ritz vectors corresponding to K harmonic Ritz values 
%                  of smallest magnitude
%         H2       upper Hessenburg matrix computed GCRODR relations
%         V        N-by-M+1 matrix containing Krylov basis computed by GMRES
%         U        basis for recycled subspace
%         C        C = A*U, as per GCRODR relations
% OUTPUT: HARMVECS basis for span of K harmonic Ritz vectors
function harmVecs =lgetHarmVecs2(m,k,H2,V,U,C)
B = H2' * H2;

% A = | C'*U        0 |
%     | V_{m+1}'*U  I |
A = zeros(m+1,m);
A(1:k,1:k) = C' * U;
A(k+1:m+1,1:k) = V' * U;
A(k+1:m,k+1:m) = eye(m-k);
A = H2' * A;

N=zeros(m);
N(1:k,1:k)=U'*U;
N(k+1:m,1:k)=A(k+1:m,1:k);
N(1:k,k+1:m)=A(k+1:m,1:k)';
N(k+1:m,k+1:m)=eye(m-k);

% Compute k smallest refined harmonic Ritz pairs.
[vec,val]= eig(B,A);
ritz_1=diag(val);
[ritzz,ir]=sort(abs(ritz_1));
for i=1:k
  ritz(i)=ritz_1(ir(i));
  W_m=B+abs(ritz(i))^2*N-(ritz(i)*A+ritz(i)'*A');
 [v ld u]=svd(W_m);
  harmVecs(:,i)=u(:,m);
end


  