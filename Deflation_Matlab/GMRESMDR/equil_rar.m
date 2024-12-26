% April 2011
% by QU Long

function [eqA,R,C] = equil_rar (A)
% A = rand(5) ;
thresh = 0.1;
[m, n] = size(A);

R = max(abs(A'));    %Find the maximum element in each row
rcmax = max(R);
rcmin = min(R);

if rcmin == 0.
    fprintf('The rank of the matrix is not sufficient!\n');
    fprintf('System is not inversible!\n');
    return;
end

R = sqrt(1 ./ R);    % Invert the scale factors
eqA = A;

if (rcmin/rcmax) < thresh
    fprintf('RAR scaling : yes!\n');
    eqA = spdiags(R',0,m,m)*eqA*spdiags(R',0,m,m); %1./sqrt(R)*A*1./sqrt(R) . 1/sqrt(R_i)*a_ij*1/sqrt(R_j)
else
    fprintf('RAR scaling : no\n');
    R = ones(1,m);
end

C = R;
return;

