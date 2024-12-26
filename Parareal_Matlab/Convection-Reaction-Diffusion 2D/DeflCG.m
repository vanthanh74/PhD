% This program test the effect of higher order projection vectors
% and the deflation approach
%
% programmer: Kees Vuik
% e-mail    : c.vuik@math.tudelft.nl
% date      : 06-10-2005

function [U1] = DeflCG(A,b,U0)

% b = U0 + F;


n = size(A,1); % dimension of matrix
m = sqrt(n);   % number of blocks

% define the cubic projection vectors
   for i = 1: m
      Z(n,i+3*m) = 0;
      for j = 1:n/m
         Z(j+(i-1)*n/m,i+3*m) = j*j*j;
      end
   end


Z  = sparse(Z);
AZ = sparse(A*Z);
E  = Z'*AZ;
P = sparse(eye(n)-AZ*inv(E)*Z');

% deflated cg

xtrue = A\b;
acc = 1e-14;
u = Z*inv(E)*Z'*b;
normxtrue = norm(xtrue);

i = 1;
U0 = 0*b;
r = P*b;
p = r;

residu(i) = norm(r);
fout(i)   = norm(xtrue-u)/normxtrue;

while  (i < 400) && (residu(i) > acc)
   i = i+1;
   PAp = P*(A*p);
   alpha = (r'*r)/(p'*PAp);
   U1 = U0+alpha*p;
   U0 = U1;  
   r = r-alpha*PAp;
   beta = (r'*r)/(residu(i-1)^2);
   p = r+beta*p;
   residu(i) = norm(r);
   fout(i) = norm(xtrue-(u+P'*U1))/normxtrue;
end

