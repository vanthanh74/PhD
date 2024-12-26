%   (I-Q2Q2t + Q2*R22^-1*Q2')*A_kND*x = (I-Q2Q2t + Q2*R22^-1*Q2')*b
function y = defpreconditioner(Q2,R22,x)
% aa = A*x;
% bb = Q2'*aa;
bb = Q2'*x;
cc = R22\bb;
% cc = pinv(R22)*bb;
dd = Q2*cc;
ee = Q2*bb;
y  = x - ee + dd;
% y = y*x;
    %y=diag(Pr)\x;
end

% function y = defpreconditioner(Q2,R22,maxSV,x)
% % aa = A*x;
% % bb = Q2'*aa;
% bb = Q2'*x;
% cc = R22\bb;
% % cc = pinv(R22)*bb;
% dd = maxSV*Q2*cc;
% ee = Q2*bb;
% y  = x - ee + dd;
% % y = y*x;
%     %y=diag(Pr)\x;
% end