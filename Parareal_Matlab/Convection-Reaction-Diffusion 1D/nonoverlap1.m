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
C = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0;];
R1 = [1 0 0 0 0];
% D1 = [1 0 0;0 1/2 0;0 0 1];
A1 = R1*A*R1';

R2 = [0 1 0 0 0;0 0 1 0 0];
A2 = R2*A*R2';

R3 = [0 0 0 1 0;0 0 0 0 1];
A3 = R3*A*R3';

Rf = [R1;R2;R3];
Rf*A*Rf';


A_bar = R1'*A1*R1 + R2'*A2*R2 + R3'*A3*R3;
A_inverse_bar = R1'*A1^-1*R1 + R2'*A2^-1*R2 + R3'*A3^-1*R3;
A_mul = A_bar\A;
A_bar*A_mul;
B_mul = [0 0 0 0 0;0 0 0 0 0; phi^2 -phi 0 0 0;0 0 0 0 0;0 0 phi^2 -phi 0];
C_mul = [0 0 0 0 0; 0 0 0 0 0;-phi^2 0 0 0 0; 0 0 0 0 0; 0  0 -phi^2 0 0];
D_mul = [0 0 0 0 0; 0 0 0 0 0;0 -phi  0 0 0; 0 0 0 0 0; 0 0 0 -phi 0];
A-B_mul;
C_mul+D_mul;
% additive
M1  = ( R1'*A1*R1 + R2'*A2*R2 + R3'*A3*R3 ) + (R'*B*R - R'*R);
M1^-1;

% multiplicative-like
M2  = ( R1'*A1*R1 + R2'*A2*R2 + R3'*A3*R3 )*(R'*B*R + eye(n)-R'*R)
M2_inverse =M2^-1

% multiplicative
M3 =  (R1'*A1*R1 + R2'*A2*R2 + R3'*A3*R3) + (R'*B*R + eye(n)-R'*R ) ...
        - (R1'*A1*R1 + R2'*A2*R2 + R3'*A3*R3)*A*(R'*B*R  + eye(n)-R'*R )
        
M3_inverse =  (R1'*A1^-1*R1 + R2'*A2^-1*R2 + R3'*A3^-1*R3) + (R'*B^-1*R + eye(n)-R'*R) ...
        - (R'*B^-1*R+ eye(n)-R'*R )*A*(R1'*A1^-1*R1 + R2'*A2^-1*R2 + R3'*A3^-1*R3)
             
(R1'*A1*R1 + R2'*A2*R2 + R3'*A3*R3)*(R'*B*R);
% (R'*B^-1*R + eye(n)-R'*R)*A*( R1'*A1^-1*R1 + R2'*A2^-1*R2 + R3'*A3^-1*R3)*A

% ( -R'*R)*A*( R1'*A1^-1*R1 + R2'*A2^-1*R2 + R3'*A3^-1*R3 - R'*R)
P = (R1'*A1^-1*R1 + R2'*A2^-1*R2 + R3'*A3^-1*R3);
C = R'*B*R + eye(n)-R'*R ;
C1 = R'*B^-1*R + eye(n)-R'*R ;
R*C*R';
C^-1;
R'*B^-1*R;
% % additive
% M2  = R'*B1*R + ( R1'*A1*R1 + R2'*A2*R2 + R3'*A3*R3 )
% M2^-1
% % R'*B1*R + eye(n)
% 
% (R'*B*R + eye(n) - R'*R)^-1
% % (R'*(B^-1)*R + eye(n)) - (R'*(B^-1)*R*R'*R*R'*(B^-1)*R)/(1+R*R'*(B^-1)*R*R')
% A_tilde = R'*B*R + eye(n);
% % sherman-morrison formula
% A_tilde^-1 - A_tilde^-1*(-R')* ( eye(3) + R*A_tilde^-1*(-R') )^-1 *R*A_tilde^-1

%F-relaxation
F = eye(5) - C1*P*A
