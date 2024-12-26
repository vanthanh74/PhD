clear
clc
close all


n = 100;
k = 20;
kappa = 1;
[A,b,x] = heat(n,kappa);
% A= rand(n,n);

% Matlab QR with Column Pivoting
[Q,R,P] = qr(A,0);

% R22_ref= Q(:,P(end-k+1:end))*R(P(end-k+1:end),:);
R22_ref = R(end-k+1:end,end-k+1:end);
[U_ref,S_ref,V_ref] = svd(R22_ref);
[U,S,V] = svd(A);

S_diag=diag(S);
S_ref;
S_A=S_diag(end-k+1:end);
norm(S_ref-S_A,2);

S_A_last=S_A;
% QR with tournament pivoting
% A1= zeros(k,k);
% A2= zeros(k,k);

A1= A(:,1:n/2);
A2= A(:,n/2+1:end);

%QRCP for A1
[Q1,R1,P1] = qr(A1,0);

%Select the last k columns from A1
P1=sort(P1(end-k+1:end));
A11 =  A1(:,P1);

%QRCP for A2
[Q2,R2,P2] = qr(A2,0);

%Select the last k columns from A2
P2=sort(P2(end-k+1:end));
A22 =  A2(:,P2);

%Put the last k columns from A1 and A2 together
% P_TP = sort([P1(end-k/2+1:end) P2(end-k/2+1:end)])
A3 = [A11 A22];

[Q3,R3,P3] = qr(A3,0);
P31=sort(P3(end-k+1:end));
P30=sort(P3(1:end-k));
A33 =  A3(:,P31);

% Put the worst columns of A to the end
A_permuted= [A(:,P30) A(:,P31)]

[Q4,R4,P4] = qr(A_permuted,0);

R22_TP = R4(end-k+1:end,end-k+1:end);


% R22_PT = A3(end-k+1:end,end-k+1:end)

[U_TP,S_TP,V_TP] = svd(R22_TP);
S_TP_dia=diag(S_TP);
S_TP_last = S_TP_dia(end-k+1:end);
S_A;
S_A_last
S_TP_last
norm(S_A_last-S_TP_last,2)

figure
semilogy(1:k,S_A_last,'r',1:k,S_TP_last,'b')
legend('Singular Values of A','Singular Values of A_3')