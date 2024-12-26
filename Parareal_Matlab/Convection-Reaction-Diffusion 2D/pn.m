function p = pn(n,tau)

p = 1;
for i = 1:n   
    p = p*(1-tau^i);
end

p = 1/p;
