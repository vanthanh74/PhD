clear all
close all
clc

matrix='CurlCurl_0.mtx';  %sherman3 N=5005 nz=20033,   saylr4
A =  mmread(matrix);
fprintf('Size of A - %d\n', length(A));
fprintf('Non zeros of A - %d\n', nnz(A));

n = length(A);
npes =    8;  %number of subdomains nested disseciton npes=2^n
k= 20;   % number of smallest singular values to deflate
% maxit = 25;
tol = 1e-6;

%A_ND = A(p,p),     [p,ip,sizes] = metismex('NodeNDP',A,npes); %
% matrix A_kND obtain from QRCP with nested dissection with the last k
% columns corresponding to the k smallest singular values
[A_ND,A_kND,P_kND]= NDQRCPSingularValuesApproximationPk(A,npes,k);


% A = randi([1 8],8)
% I=eye(size(A));
% p=[1 3 5 7 2 4 6 8];
% I(:,p)
% A*I(:,p)

I=eye(n,n);
P_kND_mat = I(:,P_kND);
% norm(A_ND*P_kND_mat-A_kND)
% norm(A_ND-A_kND*P_kND_mat')

% A_kND*P_C = Q_final*R_final=  [Q1 Q2] *   [R11  R12]
%                                           [     R22]
% A_kND = full(A_kND);
[Q_final,R_final] = qr(full(A_kND),'matrix');   %   [Q,R,E] = QR(A) produces unitary Q, upper triangular R and a
%   permutation matrix E so that A*E = Q*R. The column permutation E is
%   chosen to reduce fill-in in R.
%----------------------------------------------------
% A_permutted = A_final;
P_C = P_kND_mat;
% normAPC_QR=normest(A_kND*P_C-Q_final*R_final)
Q1=Q_final(1:n,1:n-k);  % size n x n-k
Q2=Q_final(1:n,end-k+1:end); % size n x k
R11=R_final(1:n-k,1:n-k); % size n-k x n-k
R12=R_final(1:n-k,n-k+1:end); % size n-k x k
R22=R_final(end-k+1:end,end-k+1:end); % size k x k
normR22=norm(R22)
%%Approximate for the right null space of A
V2_tilde = P_C*[-R11\R12;eye(k)];  % size n x k
Q2Q2t = Q2*Q2';
Q2tQ1 = Q2'*Q1;
PctPc=P_C'*P_C;
PcPct=P_C*P_C';
% zerosIdentity=R22^-1*Q2'*A_kND;   % [0 I] size k x n
% aaa=V2_tilde*V2_tilde^-1;
%verify
% Q22 = A_kND*V2_tilde*R22^-1;
% V2tilde=A_kND^-1*(Q2*R22);
% AV2tilde=A_kND*V2_tilde;
% Q2R22=Q2*R22;
% norm(AV2tilde-Q2R22)



% %------------------------------------------------------------------
% M^-1 * A_kND
% (I-Q2Q2t)*A_kND*x = (I-Q2Q2t)*b

MinvA = (I-Q2*Q2'+V2_tilde*(Q2'*A_ND*V2_tilde)^-1*Q2')*A_kND;
rhsMinvA=[Q1 V2_tilde]*[R11 R12; zeros(k,n-k) eye(k,k)];%*P_C'
normMinvA=norm(MinvA-rhsMinvA)
sA=svd(full(A_kND));
s=svd(MinvA);
sA_k = sA(end-k:end);
s_k = s(end-k:end);

% 
% aa=R22^-1*Q2'*A_kND - zerosIdentity*P_C';
% aa=V2_tilde*zerosIdentity*P_C';
% bb=aa*P_C';
% % ------------------------------------------------------------------
% M^-1 * A_kND  
% (I-Q2Q2t + Q2*R22^-1*Q2')*A_kND*x = (I-Q2Q2t + Q2*R22^-1*Q2')*b

MinvA1 = (I-Q2*Q2' +  Q2*R22^-1*Q2')*A_kND;
rhsMinvA1=[Q1 Q2]*[R11 R12; zeros(k,n-k) eye(k,k)]; %*P_C'
normMinvA1=normest(MinvA1-rhsMinvA1)
s1=svd(MinvA1);
s_k1 = s1(end-k:end);

%%
% right-hand side
b=rand(n,1);
b=b/norm(b);
x_exact = A_kND\b;
%%
%scaling 
[A_kND_scaling,R00,C00] = equil_rar(A_kND);


%--------------------------------------------------------
%   (I-Q2Q2t + Q2*R22^-1*Q2')*A_kND*x = (I-Q2Q2t + Q2*R22^-1*Q2')*b
aa = A_kND_scaling;
bb = Q2'*aa;
cc = R22\bb;
dd = Q2*cc;
ee = Q2*bb;
y  = aa - ee + dd;

%---------------------------------------------------------
% (I-Q2Q2t)*A_kND*x = (I-Q2Q2t)*b
% bb = Q2'*A_kND_scaling;
% cc = Q2*bb;
% y = A_kND_scaling - cc;
%---------------------------------------------------------
M  = @(x) y*x;%preconditioner; 

x = gmres(A_kND_scaling,b,[],1e-6,2600,M);  %