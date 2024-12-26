function U = mdl1(u0,T,n,type)
% backward Euler
U=[u0];
dt = T/n;

if strcmp(type,'coarse')
    if n == 1
       U(2) =  1/(2*dt) - sqrt( 1/(4*dt^2) - U(1)/dt );
    else
        for i = 2:n+1
             U(i) =  1/(2*dt) - sqrt( 1/(4*dt^2) - U(i-1)/dt );
        end
    end
elseif strcmp(type,'fine')
    for i = 2:n+1
         U(i) =  1/(2*dt) - sqrt( 1/(4*dt^2) - U(i-1)/dt );
    end
end