%   (I-Q2Q2t + Q2*R22^-1*Q2')*A_kND*x = (I-Q2Q2t + Q2*R22^-1*Q2')*b
%   with scaling A_kND
aa = A_kND_scaling;
bb = Q2'*aa;
cc = R22\bb;
dd = Q2*cc;
ee = Q2*bb;
y  = aa - ee + dd;
%---------------------------------------------------------------------------------
M  = @(x) y*x;%preconditioner; 
x = gmres(A_kND_scaling,b,[],1e-6,2600,M);  
%---------------------------------------------------------------------------------
% saylr4 , size(A_kND) = 3564x3564 RAR scaling : yes!
gmres stopped at iteration 2600 without converging to the desired tolerance 1e-06
because the maximum number of iterations was reached.
The iterate returned (number 2600) has relative residual 0.041.
%---------------------------------------------------------------------------------
% sherman3, size(A_kND) = 5005x5005 RAR scaling : yes!
gmres converged at iteration 1621 to a solution with relative residual 9.9e-07.
%---------------------------------------------------------------------------------
% lshp3466, size(A_kND) = 3466x3466 RAR scaling : no
gmres stopped at iteration 2600 without converging to the desired tolerance 1e-06
because the maximum number of iterations was reached.
The iterate returned (number 2600) has relative residual 0.035.
%---------------------------------------------------------------------------------
% CurlCurl_0, size(A_kND) = 11083x11083 RAR scaling : yes!
gmres converged at iteration 53 to a solution with relative residual 9.3e-07.
%---------------------------------------------------------------------------------
% bcsstm25, size(A_kND) = 15439x15439 RAR scaling : yes!
gmres stopped at iteration 2600 without converging to the desired tolerance 1e-06
because the maximum number of iterations was reached.
The iterate returned (number 2600) has relative residual 0.16.