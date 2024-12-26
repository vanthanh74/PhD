function [U1] = BackwardEulerCoarse_2DConvDiff(kappa,dX_coarse,dY_coarse,dT,f,u0,t0,nx_coarse,ny_coarse)

global xC yC p1 p2 p3


% u0 = u0_restriction
% t0 = 0

U0 = zeros((nx_coarse-1)*(ny_coarse-1),1);
for i = 1:nx_coarse-1
    for j=1:ny_coarse-1
        U0((i-1)*(nx_coarse-1)+j) = u0(i+1,j+1);
    end    
end

F  = zeros((nx_coarse-1)*(ny_coarse-1),1);
for i=1:nx_coarse-1
    for j=1:ny_coarse-1
        F((i-1)*(nx_coarse-1) +j)= dT*f(i*dX_coarse,j*dX_coarse,t0+dT);
    end
end


r1 = kappa*dT/(dX_coarse^2);
k1 = p1*dT/dX_coarse;

r2 = kappa*dT/(dY_coarse^2);
k2 = p2*dT/dY_coarse;

A = -(r1+k1)*diag(ones((ny_coarse-1)*(ny_coarse-1) - 1,1),-1) ...
    + (1+2*r1+2*r2 - p3*dT)*eye((ny_coarse-1)*(ny_coarse-1)) ...
    - (r1-k1)*diag(ones((ny_coarse-1)*(ny_coarse-1)-1,1),1) ...
    - (r2+k2)*diag(ones((ny_coarse-1)*(ny_coarse-1)-(ny_coarse-1),1),-(ny_coarse-1)) ...
    - (r2-k2)*diag(ones((ny_coarse-1)*(ny_coarse-1)-(ny_coarse-1),1),(ny_coarse-1)) ;

B = zeros((nx_coarse-1)*(ny_coarse-1),(nx_coarse-1)*(ny_coarse-1));
    for i = 1:nx_coarse-2
        B(i*(nx_coarse-1),i*(nx_coarse-1)+1) = r2-k2;
        B(i*(nx_coarse-1)+1,i*(nx_coarse-1)) = r2+k2;
    end
    
A = A+B;
A = sparse(A);

% I = r2*eye(ny_coarse-1,ny_coarse-1);

% A = zeros((nx_coarse-1)*(ny_coarse-1),(nx_coarse-1)*(ny_coarse-1))

% U1 = A\(U0 + F);
   U1 = cgs(A,U0+F,1e-8,50);

U_matrix = zeros((nx_coarse+1),(ny_coarse+1));
for i = 2:nx_coarse
    for j =2:ny_coarse
         U_matrix(i,j) =  U1((j-2)*(nx_coarse-1)+i-1) ;
    end
end
U1 = U_matrix;
