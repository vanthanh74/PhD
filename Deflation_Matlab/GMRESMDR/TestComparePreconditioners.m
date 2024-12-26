clear all
close all
clc

matrix='saylr4.mtx';  %sherman3 N=5005 nz=20033,   saylr4
A =  mmread(matrix);
fprintf('Size of A - %d\n', length(A));
fprintf('Non zeros of A - %d\n', nnz(A));


%scaling 
[A,R00,C00] = equil_rar(A);

n = length(A);
npes =    8;  %number of subdomains nested disseciton npes=2^n
k=20;   % number of smallest singular values to deflate
% maxit = 25;
tol = 1e-6;

%A_ND = A(p,p),     [p,ip,sizes] = metismex('NodeNDP',A,npes); %
% matrix A_kND obtain from QRCP with nested dissection with the last k
% columns corresponding to the k smallest singular values
[A_ND,A_kND,P_kND,p]= NDQRCPSingularValuesApproximationPk(A,npes,k);


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
% V2_tilde = P_C*[-R11^-1*R12;eye(k)];  % size n x k
% Q2Q2t = Q2*Q2';
% Q2tQ1 = Q2'*Q1;
% PctPc=P_C'*P_C;
% PcPct=P_C*P_C';
% zerosIdentity=R22^-1*Q2'*A_kND;   % [0 I] size k x n
% aaa=V2_tilde*V2_tilde^-1;
%verify
% Q22 = A_kND*V2_tilde*R22^-1;
% V2tilde=A_kND^-1*(Q2*R22);
% AV2tilde=A_kND*V2_tilde;
% Q2R22=Q2*R22;
% norm(AV2tilde-Q2R22)



%------------------------------------------------------------------
% % M^-1 * A_kND
% MinvA = (I-Q2Q2t+V2_tilde*(Q2'*A_kND*V2_tilde)^-1*Q2')*A_kND;
% rhsMinvA=[Q1 V2_tilde]*[R11 R12; zeros(k,n-k) eye(k,k)]*P_C';
% normMinvA=norm(MinvA-rhsMinvA);
% sA=svd(full(A_kND));
% s=svd(MinvA);
% sA_k = sA(end-k:end);
% s_k = s(end-k:end);


% aa=R22^-1*Q2'*A_kND - zerosIdentity*P_C';
% aa=V2_tilde*zerosIdentity*P_C';
% bb=aa*P_C';
%------------------------------------------------------------------
% % M^-1 * A_kND
% MinvA1 = (I-Q2Q2t +  Q2*R22^-1*Q2')*A_kND;
% rhsMinvA1=[Q1 Q2]*[R11 R12; zeros(k,n-k) eye(k,k)]*P_C';
% normMinvA1=norm(MinvA1-rhsMinvA1);
% s1=svd(MinvA1);
% s_k1 = s1(end-k:end);


%We know that $V2_tilde$ is close to $V_2$ in the sense that $V2_tilde$ represents well the smallest singular values
% [U,Sigma,V] = svd(full(A));

% [Q,R,P_C] = qr(full(A),'matrix');

% V2= V(:,end-k+1,end); % right singular vectors associated to the k smallest singular values
% normV_V2tildenormest=normest(V(:,end-k+1,end)-V2_tilde)
% size(V)
% size(V2_tilde)

% right-hand side
b=rand(n,1);
b=b/norm(b);

x_exact = A\b;



nb_subdomains = 8;
% % Partitioning using Metis
[p1, edgecut, s, begin_in, end_in] = DefineAlpha(A, nb_subdomains);
A = A(p1, p1);
for i = 1 : nb_subdomains
    Pr{i} = A(begin_in(i) : end_in(i), begin_in(i) : end_in(i));
end


%--------------------------------------------------------
%   (I-Q2Q2t + Q2*R22^-1*Q2')*A_kND*x = (I-Q2Q2t + Q2*R22^-1*Q2')*b

M1 = @(x) defpreconditioner(Q2,R22,x);   % deflation preconditioner
M2 = @(x) blcLS(Pr, begin_in, end_in, x); % block jacobi preconditioner
M3  = @(x)  defpreconditioner(Q2,R22,x) + blcLS(Pr, begin_in, end_in,x);  
%---------------------------------------------------------
maxit=500;
tol=1e-6;
%gmres with deflation preconditioner
[x_gmres1,flag1,relres1,iter1,resvec1]=gmres(A,b,[],tol,maxit,M1);
iter1
relres1=relres1
%--------------------------------------------------------
%gmres with deflation preconditioner
[x_gmres2,flag2,relres2,iter2,resvec2]=gmres(A,b,[],tol,maxit,M2);
iter2
relres2=relres2
%--------------------------------------------------------
%gmres with mixed preconditioner
[x_gmres3,flag3,relres3,iter3,resvec3]=gmres(A,b,[],tol,maxit,M3);
iter3
relres3=relres3
 
% norm(x_gmres2-x_exact)
% figure
% plot(1:n,x_gmres2,'db',1:n,x_exact,'.r')



