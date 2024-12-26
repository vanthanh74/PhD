
% Ex1-RST

function  [ A b x ]  = ex2_rst( m )

m1 = 2^round(log2(m));
m2 = 2^round(log2(m/12));
m3 = 2^round(log2(m/20));
mmm = [m1 12*m2 20*m3];
[m ind] = min(abs(mmm - m));
m = mmm(ind);
n = m/2;

K = 10;
s_ = 1e-6;

j = 1 : K;
sigma(j) = s_ .^(floor(j/2)/5);

j = (K+1) : n;
sigma(j) = s_ * (n-j)./(n-K);

U = hadamard( m )/sqrt(m);
V = hadamard( n )/sqrt(n);

A = U * [diag(sigma); zeros(m-n,n)] * V';

x = ones(n,1);
b = A * x;
