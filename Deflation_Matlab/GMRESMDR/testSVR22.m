clear all
close all
clc

matrix='saylr4.mtx';  %sherman3 N=5005 nz=20033,  bcsstm25 bcsstk35  saylr4 lshp3466 CurlCurl_0 crystk03.mtx ct20stif.mtx bcsstm39.mtx
A =  mmread(matrix);

% n = 15000;
% A = sprandn(n,n,0.00001);
% spy(A)
fprintf('Size of A - %d\n', length(A));
fprintf('Non zeros of A - %d\n', nnz(A));
% spy(A)
% svdA=svd(full(A));
%scaling 
[A,R00,C00] = equil_rar(A);
% svdAscaling=svd(full(A));
% A0=A;

n = length(A);
npes =    8;  %number of subdomains nested disseciton npes=2^n
k=20;   % number of smallest singular values to deflate
% maxit = 25;
tol = 1e-6;

%A_ND = A(p,p),     [p,ip,sizes] = metismex('NodeNDP',A,npes); %
% matrix A_kND obtain from QRCP with nested dissection with the last k
% columns corresponding to the k smallest singular values
[A_ND,A_kND,P_kND,p]= NDQRCPSingularValuesApproximationPk(A,npes,k);

%----------------------------------------------------
%verify
% A = randi([1 8],8)
% I=eye(size(A));
% p=[1 3 5 7 2 4 6 8];
% I(:,p)
% A*I(:,p)

I=speye(n,n);
P_ND = I(:,p); %P_ND
Pi = I(:,P_kND); % Pi
% norm(A_ND*P_kND_mat-A_kND)
% norm(A_ND-A_kND*P_kND_mat')
% normest(A_kND-P_ND'*A*P_ND*Pi)
% A_kND*P_C = Q_final*R_final=  [Q1 Q2] *   [R11  R12]
%                                           [     R22]


%----------------------------------------------------
% k-way partition
nb_subdomains = 8;
% % Partitioning using Metis
[p1, edgecut, s, begin_in, end_in] = DefineAlpha(A, nb_subdomains);
P_kw = I(:,p1);    %  P_kw
A = A(p1, p1);
for i = 1 : nb_subdomains
    Pr{i} = A(begin_in(i) : end_in(i), begin_in(i) : end_in(i));
end
%----------------------------------------------------
% figure
% spy(A)

% A_kND = full(A_kND);
[Q_final,R_final] = qr(full(A_kND),'matrix');   
% normest(A_kND-Q_final*R_final)

Q_tilde = P_kw'*P_ND*Q_final;
% R_tilde = R_final*Pi'*P_ND'*P_kw;
% Q_tilde1 = P_kw'*P_ND'*Q_final;

% A_permutted = A_final;
P_C = Pi;
% normAPC_QR=normest(A_kND*P_C-Q_final*R_final)
Q1=Q_tilde(1:n,1:n-k);  % size n x n-k
Q2=Q_tilde(1:n,end-k+1:end); % size n x k
R11=R_final(1:n-k,1:n-k); % size n-k x n-k
R12=R_final(1:n-k,n-k+1:end); % size n-k x k
R22=R_final(end-k+1:end,end-k+1:end); % size k x k
normR22=normest(R22)

% figure
% subplot(121)
% spy(R_final)
% subplot(122)
% spy(R22)
% figure
% spy(R22)


% AA = rand(1000);
% BB = AA(end-k+1:end,end-k+1:end);
% 
% svAA = svd(AA);
% svAA_k=svAA(end-k+1:end);
% 
% svBB = svd(BB);
% svBB_k=svBB(end-k+1:end);
% 
% [svAA_k svBB_k]

%svd(A)
sA=svd(full(A));
% [u,sA,v]=svd(full(A));
k_smallest_sv_svd = sA(end-k+1:end);

% sA_kND = svd(full(A_kND));
% sA_kND_k =sA_kND(end-k+1:end);

%svd(R_final)
sR=svd(R_final);
sR_k = sR(end-k+1:end);
 
sR22=svd(R22);
sR22_k = sR22(end-k+1:end);

[k_smallest_sv_svd  sR_k sR22_k]