function result = ff(x,t)
global  a b c L nx_fine xF tF
% x=xF(2:nx_fine);
% t = tF(2);


result = b.*(exp(-2.*t).*(L - x).^2 - x.*exp(-2.*t).*(2.*L - 2.*x)) - a.*(2.*x.*exp(-2.*t) - 2.*exp(-2.*t).*(2.*L - 2.*x)) - 2.*x.*exp(-2.*t).*(L - x).^2 - c.*x.*exp(-2.*t).*(L - x).^2;
