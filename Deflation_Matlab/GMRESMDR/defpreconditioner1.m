%   (I-Q2Q2t + Q2*R22^-1*Q2')*P_ND'*x = (I-Q2Q2t + Q2*R22^-1*Q2')*b
function y = defpreconditioner1(Q2,R22,P_ND,x)
% aa = A*x;
% bb = Q2'*aa;
bb = Q2'*P_ND'*x;
cc = R22\bb;
% cc = pinv(R22)*bb;
dd = Q2*cc;
ee = Q2*bb;
y  = P_ND'*x - ee + dd;
% y = y*x;
    %y=diag(Pr)\x;
end