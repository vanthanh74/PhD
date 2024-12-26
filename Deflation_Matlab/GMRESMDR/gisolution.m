function [x1,g,g1]=gisolution(A,b)
% A=lhs;
% b=rhs;
[m,n] = size(A);
beta = 1.2 ./ max(eig (A^2*A'^2));
z0 = beta*(A*(A')^2);
x0 = A*A(:,1);
I = eye(m);
iter = 0;
f = 1;
maxiter = 20;
tic ;
k=0;
while (f > eps) && (iter <maxiter)
    k=k+1
    z1 = z0 + z0*(I - A*z0);
    x1 = x0 + z1*(b - A*x0);
    iter = iter + 1;
    g = norm(A*x0 - b, 2);
    g1 = norm(x1 - x0, 2);
    B = [g,g1];
    f = max(B);
    z0 = z1;
    x0 = x1;
end
toc;