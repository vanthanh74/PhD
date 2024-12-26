clear all
close all
clc




syms phi phiDT

R = [1 0 0 0 0;0 0 1 0 0;0 0 0 0 1];
% D = [1/2 0 0; -2/3*phiDT 1/3 0; -phiDT^2/2 -1/2*phiDT 1/2];
% D = [2 0 0 0 0; 0 1 0 0 0;0 0 3 0 0 ; 0 0 0 1 0;0 0  0 0 2];
B = [0 0 0;
    -phiDT 0 0;
    0 -phiDT 0];


A =[1 0 0 0 0 ; -phi 1 0 0 0; 0 -phi 1 0 0;
    0 0 -phi 1 0; 0 0 0 -phi 1];
n=length(A);

R1 = [1 0 0 0 0];
% D1 = [1 0 0;0 1/2 0;0 0 1];
A1 = R1*A*R1';

R2 = [0 1 0 0 0;0 0 1 0 0];
A2 = R2*A*R2';

R3 = [0 0 0 1 0;0 0 0 0 1];
A3 = R3*A*R3';

M = R'*B*R + R1'*A1*R1 + R2'*A2*R2 + R3'*A3*R3

Id = eye(n);
M1 = (Id - (Id - R1'*A1^-1*R1*A)*(Id - R2'*A2^-1*R2*A)*(Id - R3'*A3^-1*R3*A))*(Id-R'*B*R)