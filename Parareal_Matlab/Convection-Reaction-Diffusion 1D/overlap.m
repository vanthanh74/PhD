clear all
close all
clc

syms phi phiDT

R = [1 0 0 0 0;0 0 1 0 0;0 0 0 0 1];
% D = [1/2 0 0; -2/3*phiDT 1/3 0; -phiDT^2/2 -1/2*phiDT 1/2];
% D = [2 0 0 0 0; 0 1 0 0 0;0 0 3 0 0 ; 0 0 0 1 0;0 0  0 0 2];
B = [1 0 0;
    -phiDT 1 0;
    0 -phiDT 1];
R'*B*R;

A =[1 0 0 0 0 ; -phi 1 0 0 0; 0 -phi 1 0 0;
    0 0 -phi 1 0; 0 0 0 -phi 1];
n  = length(A);

R1 = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0];
D1 = [1 0 0; 0 1 0; 0 0 1/2];
A1 =R1*A*R1';

R2 = [0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1];
D2 = [1/2 0 0;0 1 0; 0 0 1];
A2 =R2*A*R2';

sumRtDR =  R1'*D1*R1 + R2'*D2*R2

% additive preconditioner
M1 = (R'*B*R -R'*R) + (R1'*D1*A1*R1 + R2'*D2*A2*R2)

M2 = (R'*B*R + eye(n) -2*R'*R) + (R1'*D1*A1*R1 + R2'*D2*A2*R2)


% multiplicative preconditioner
M3 = (R1'*D1*A1*R1 + R2'*D2*A2*R2)*(R'*B*R + eye(n)-R'*R) 

M4 = (R1'*D1*A1*R1 + R2'*D2*A2*R2)*(R'*B*R + 2*(eye(n)-R'*R) )
% M4 = (R'*B*R + eye(n)-R'*R )*(R1'*D1*A1*R1 * R2'*D2*A2*R2)



%-----------------------------------------------------------------
% DD  = [1/2 0 0; 0 1/2 0; 0 0 1/2];
% BB = [1/2 0 0;
%     -phiDT 1/2 0;
%     0 -phiDT 1/2];
% 
% R'*BB*R
% DD1 = [1/2 0 0; 0 1 0; 0 0 1/4];
% DD2 = [1/4 0 0;0 1 0; 0 0 1/2];
% 
% sumRtDR = R'*DD*R  + R1'*DD1*R1 + R2'*DD2*R2
% A11 = R1*A*R1';
% A22 = R2*A*R2';
% M = R'*B*R + R1'*DD1*A11*R1 + R2'*DD2*A22*R2