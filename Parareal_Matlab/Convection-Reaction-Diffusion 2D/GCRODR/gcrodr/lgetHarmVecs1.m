% getHarmVecs1      For use with GCRODR
%
% Determines harmonic Ritz vectors using matrix H computed from a
% GMRES iteration. For this case, the harmonic Ritz values are the
% eigenvalues of H
% 
% INPUT:  M        dimension of upper Hessenburg matrix H
%         K        select and return basis for space spanned by K harmonic 
%                  Ritz vectors corresponding to K harmonic Ritz values 
%                  of smallest magnitude
%         H        M+1-by-M upper Hessenburg matrix computed from GMRES 
% OUTPUT: HARMVECS basis for span of K harmonic Ritz vectors

function harmVecs = lgetHarmVecs1(m,k,H)

% Build matrix for eigenvalue problem.
harmRitzMat = H(1:m,:)' \ speye(m);
harmRitzMat(1:m,1:m-1) = 0;
harmRitzMat = H(1:m,:) + H(m+1,m)^2 * harmRitzMat;
H_m=H'*H;
% Compute k smallest harmonic Ritz pairs.
[vec, val] = eig(H_m,H(1:m,:)');
ritz_1=diag(val);
% Sort by magnitide of eigenvalue
[ritzz,ir]=sort(abs(ritz_1));
I_m=eye(m+1,m);
for i=1:k
  ritz(i)=ritz_1(ir(i));
  W_m=H-ritz(i)*I_m;
   [v ld u]=svd(W_m);
  harmVecs(:,i)=u(:,m);
end