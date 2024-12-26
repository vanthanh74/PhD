 clear all
close all
clc


% non-overlapping subdomains

syms phi phiDT

R  = [1 0 0 0 0;0 0 1 0 0;0 0 0 0 1];
% D = [1/2 0 0; -2/3*phiDT 1/3 0; -phiDT^2/2 -1/2*phiDT 1/2];
% D = [2 0 0 0 0; 0 1 0 0 0;0 0 3 0 0 ; 0 0 0 1 0;0 0  0 0 2];
B  = [1 0 0;
    -phiDT 1 0;
    0 -phiDT 1];

B1 = [0 0 0;
    -phiDT 0 0;
    0 -phiDT 0];


A  =[1 0 0 0 0 ; -phi 1 0 0 0; 0 -phi 1 0 0;
    0 0 -phi 1 0; 0 0 0 -phi 1];
n  = length(A);

A_rearr = [1 0 -phi 0 0 ; 0 1 0 -phi 0; 0 0 1 0 0;
    -phi 0 0 1 0; 0 -phi 0 0  1];


S  = [1 0 0;
    -phi^2 1 0;
    0 -phi^2 1];

R1 = [1 0 0 0 0];
R11 = [1 0 0 0 0; 0 1 0 0 0];
% D1 = [1 0 0;0 1/2 0;0 0 1];
A1 = R1*A*R1';
A11 = R11*A*R11';

R2 = [0 1 0 0 0;0 0 1 0 0];
R22 = [0 0 1 0 0; 0 0 0 1 0];
A2 = R2*A*R2';
A22 = R22*A*R22';

R3 = [0 0 0 1 0;0 0 0 0 1];
R33 = [0 0 0 0 1];
A3 = R3*A*R3';
A33 = R33*A*R33';

Rf = [R1;R2;R3];
Rf*A*Rf'

% A_bar = R1'*A1*R1 + R2'*A2*R2 + R3'*A3*R3;
% A_mul = A_bar\A;
% A_bar*A_mul;
% B_mul = [0 0 0 0 0;0 0 0 0 0; phi^2 -phi 0 0 0;0 0 0 0 0;0 0 phi^2 -phi 0];
% C_mul = [0 0 0 0 0; 0 0 0 0 0;-phi^2 0 0 0 0; 0 0 0 0 0; 0  0 -phi^2 0 0];
% D_mul = [0 0 0 0 0; 0 0 0 0 0;0 -phi  0 0 0; 0 0 0 0 0; 0 0 0 -phi 0];
% A-B_mul;
% C_mul+D_mul;
% E=(R'*S*R + eye(n)-R'*R)\A_mul
% (R'*S*R+ eye(n)-R'*R)*E
% 
% 
% decomposition of A
% A_bar*(R'*S*R+ eye(n)-R'*R)*E

A_bar = ( R1'*A1*R1 + R2'*A2*R2 + R3'*A3*R3 )*(R'*S*R + eye(n)-R'*R)
A_bar_approx=( R1'*A1*R1 + R2'*A2*R2 + R3'*A3*R3 )*(R'*B*R + eye(n)-R'*R)

% AA = ( R1'*A1*R1 + R2'*A2*R2 + R3'*A3*R3 )*(R'*S*R + eye(n)-R'*R)*(R11'*A11*R11 + R22'*A22*R22 + R33'*A33*R33)

eye(n) - ((R'*B*R + eye(n)-R'*R)^-1)*(R'*S*R + eye(n)-R'*R)

M =  ( R1'*A1*R1 + R2'*A2*R2 + R3'*A3*R3 )*(R'*B*R + eye(n)-R'*R)
a = [phi 0 0; 0 phi 0;eye(3)]
b = a;
b(1) = 0;
b(2,2) = 0;
K=a*(eye(3) - (B^-1)*S)*b'
K\(eye(5) - (M^-1)*A)
(eye(5) - (M^-1)*A)\K