%   (I-Q2Q2t + Q2*R22^-1*Q2')*A_kND*x = (I-Q2Q2t + Q2*R22^-1*Q2')*b
function y = mixedAdditivePreconditioner(Q2,R22,Pr, begin_in, end_in,x)
% aa = A*x;
% bb = Q2'*aa;

%Deflation preconditioner
bb = Q2'*x;
cc = R22\bb;
% cc = pinv(R22)*bb;
dd = Q2*cc;
ee = Q2*bb;
y1  = x - ee + dd;

%Block Jacobi preconditioner
n = length(begin_in);
for i = 1 : n
%     if size(Pr{i},1)~=0
         y2(begin_in(i) : end_in(i), :) = Pr{i}\ x(begin_in(i) : end_in(i), :);%pinv(full(Pr{i}))* x(begin_in(i) : end_in(i), :);%
%          y(begin_in(i) : end_in(i), :) = pinv(full(Pr{i}))* x(begin_in(i) : end_in(i), :);%Pr{i}\ x(begin_in(i) : end_in(i), :);%
% 
%     end
end

% combination
y = y1 + y2;

% y = y*x;
    %y=diag(Pr)\x;
end