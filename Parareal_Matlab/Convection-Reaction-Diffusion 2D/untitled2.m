clear all
clc

tau = 0.1:0.1:1.4;
DT = 0.5;
Nt = 20;
DTn = zeros(Nt,1);
condVec = [];
for j = 1:length(tau)
    S = eye(Nt-1,Nt-1);
    for i = 2:Nt-1
       S = S + pn(i-1,tau(j))*diag(ones(Nt-i,1),-(i-1));
    end
    condVec(j,1) = cond(S)


% for n = 1:Nt
%     DTn(n) = DT*tau(j)^(n-1) 
% 
%     
%     
%     
%     
% end
    

end