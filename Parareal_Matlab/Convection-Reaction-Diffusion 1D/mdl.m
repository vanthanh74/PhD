function U = mdl(u0,T,n,type)
% backward Euler
U=[u0];
global a ;
dt = T/n;

if strcmp(type,'coarse')
    if n == 1
        U(2) = U(1)/(1-a*dt);
    else
        for i = 2:n+1
            U(i) =  U(i-1)/(1-a*dt);
        end
    end
elseif strcmp(type,'fine')
    for i = 2:n+1
        U(i) =  U(i-1)/(1-a*dt);
    end
end