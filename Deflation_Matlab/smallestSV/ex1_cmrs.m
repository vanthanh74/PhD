
% Ex1-RST

function  [ A b x ]  = ex1_cmrs( m )

n = floor(m/2);

for j = 1:m
    for k = 1:m
        C(j,k) = exp(pi*(2*j-1)/(4*n-2) * cos(pi*(2*k-1)/(2*n-1)) );
    end
end
[U S0 V0] = svd(C);

clear C

for j = 1:n
    for k = 1:n
        C(j,k) = exp(pi*(2*j-1)/(4*n-2) * cos(pi*(2*k-1)/(2*n-1)) );
    end
end
[U0 S0 V] = svd(C);

j = 1:n;
d = exp(-0.4 * (j-1) );
S = [ diag(d); zeros(m-n,n) ];


A = U * S * V';

x = randn(n,1);
%x = x/norm(x);

b = A * x;
