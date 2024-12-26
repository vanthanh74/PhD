clear all
close all
clc

matrix='sherman3.mtx';  %sherman3 N=5005 nz=20033,   saylr4 lshp3466 CurlCurl_0 crystk03.mtx ct20stif.mtx bcsstm39.mtx
A =  mmread(matrix);

% n = 15000;
% A = sprandn(n,n,0.00001);
% spy(A)
fprintf('Size of A - %d\n', length(A));
fprintf('Non zeros of A - %d\n', nnz(A));


%scaling 
[A,R00,C00] = equil_rar(A);

n = length(A);
npes =    8;  %number of subdomains nested disseciton npes=2^n
% k=20;   % number of smallest singular values to deflate
k_df= 50   % number of smallest singular values to deflatef % rank of R22
k = n-k_df; % rank of A11

% maxit = 25;
tol = 1e-6;

%A_ND = A(p,p),     [p,ip,sizes] = metismex('NodeNDP',A,npes); %
% matrix A_kND obtain from QRCP with nested dissection with the last k
% columns corresponding to the k smallest singular values
% [A_ND,A_kND,P_kND,p]= NDQRCPSingularValuesApproximationPk(A,npes,k);
[A_ND,A_kND,P_kND,p,last_k_columns_matrix,last_k_columns_matrix_temp]= NDQRCPSingularValuesApproximationPk(A,npes,k,k_df);


I=speye(n,n);
P_ND = I(:,p); %P_ND
Pi = I(:,P_kND); % Pi

figure
spy(A)

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

% M2 = @(x) blcLS(Pr, begin_in, end_in, x); % block jacobi preconditioner

% function y = blcLS(Pr, begin_in, end_in, x)
% %block linear solver  Pr^-1 * x  with k-way partition
% n = length(begin_in);
% for i = 1 : n
%          y(begin_in(i) : end_in(i), :) = Pr{i}\ x(begin_in(i) : end_in(i), :);
% end
% end


%----------------------------------------------------
% figure
% spy(Pr{2})


% A = randi([1 8],8)
% I=eye(size(A));
% p=[1 3 5 7 2 4 6 8];
% I(:,p)
% A*I(:,p)

% norm(A_ND*P_kND_mat-A_kND)
% norm(A_ND-A_kND*P_kND_mat')

% A_kND*P_C = Q_final*R_final=  [Q1 Q2] *   [R11  R12]
%                                           [     R22]
% A_kND = full(A_kND);
% [Q_final,R_final] = qr(full(A_kND),'matrix');   %   [Q,R,E] = QR(A) produces unitary Q, upper triangular R and a
%   permutation matrix E so that A*E = Q*R. The column permutation E is
%   chosen to reduce fill-in in R.

%strong RRQR for A_kND
f = 2;
[Q_final, R_final, P_final] = sRRQR(full(A_kND), f, 'rank', k);



%----------------------------------------------------

Q_tilde = P_kw'*P_ND*Q_final;


% A_permutted = A_final;
P_C = P_ND;
% normAPC_QR=normest(A_kND*P_C-Q_final*R_final)
% Q1=Q_tilde(1:n,1:n-k);  % size n x n-k
Q2=Q_tilde(1:n,end-k_df+1:end); % size n x k
% R11=R_final(1:n-k,1:n-k); % size n-k x n-k
% R12=R_final(1:n-k,n-k+1:end); % size n-k x k
R22=R_final(end-k_df+1:end,end-k_df+1:end); % size k x k
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

% % compare the k smallest singular values of
% % A
% sA=svd(full(A));
% sA_k = sA(end-k+1:end);
% sA_kND=svd(full(A_kND));
% sA_kkND = sA_kND(end-k+1:end);
% 
% %------------------------------------------------------------------
% % % M^-1 * A
% MinvA = (I-Q2*Q2' +  Q2*R22^-1*Q2')*A;
% s=svd(MinvA);
% s_k = s(end-k+1:end);
% %------------------------------------------------------------------
% % % (M_1^-1 + M_2^-1) * A
% MinvA1 = (I-Q2*Q2' +  Q2*R22^-1*Q2')*A +  blcLS(Pr, begin_in, end_in,A);
% s1=svd(MinvA1);
% s_k1 = s1(end-k+1:end);
% 
% compare_max=[sA_k s_k s_k1];

%%
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

% %A-DEF2 preconditioner
% QQ = Q2*(Q2'*A*Q2)^-1*Q2'; %Qfunc
% P = I - A*QQ;
% MinvA = P'*defpreconditioner(Q2,R22,A) + QQ;%;%(I-Q2*Q2' +  Q2*(R22^-1)*Q2')*A;%
% MinvA1 = defpreconditioner(Q2,R22,A) + blcLS(Pr, begin_in, end_in,A) - blcLS(Pr, begin_in, end_in,A)*A*defpreconditioner(Q2,R22,A);



%--------------------------------------------------------
%   (I-Q2Q2t + Q2*R22^-1*Q2')*A_kND*x = (I-Q2Q2t + Q2*R22^-1*Q2')*b
% M  = @(x) defpreconditioner(Q2,R22,x);% 
M1 = @(x) defpreconditioner(Q2,R22,x) ;%defpreconditioner(Q2,R22,x);   % deflation preconditioner
M2 = @(x) blcLS(Pr, begin_in, end_in, x); % block jacobi preconditioner
M3  = @(x) defpreconditioner(Q2,R22,x) + blcLS(Pr, begin_in, end_in,x) - blcLS(Pr, begin_in, end_in,A)*defpreconditioner(Q2,R22,x);
%---------------------------------------------------------
maxit=length(b);
tol=1e-6;

%non-preconditioned
[x_gmres0,flag0,relres0,iter0,resvec0]=gmres(A,b,[],tol,maxit);
iter0 = iter0(2)
relres0=relres0

%deflation
[x_gmres1,flag1,relres1,iter1,resvec1]=gmres(A,b,[],tol,maxit,M1);
iterDef = iter1(2)
relresDef=relres1

%block jacobi
[x_gmres2,flag2,relres2,iter2,resvec2]=gmres(A,b,[],tol,maxit,M2);
iterBlockJ= iter2(2)
relresBlockJ=relres2
 
%mixed
[x_gmres3,flag3,relres3,iter3,resvec3]=gmres(A,b,[],tol,maxit,M3);
iterMixed= iter3(2)
relresMixed=relres3
 
Compare_res_vec=[resvec0 resvec1  resvec2 resvec3]
Compare = [relres0 iter0; relresDef iterDef ; relresBlockJ iterBlockJ ;relresMixed iterMixed ]
% norm(x_gmres3-x_exact)
% figure
% plot(1:n,x_gmres3,'db',1:n,x_exact,'.r')
iter0
relres0
iterDef
relresDef
iterBlockJ
relresBlockJ
iterMixed
relresMixed

figure
semilogy(resvec0)
hold on
semilogy(resvec1)
semilogy(resvec2)
semilogy(resvec3)
legend('GMRES','Strong RRQR deflated GMRES','Block Jacobi preconditioned GMRES','Block Jacobi and Deflation preconditioned GMRES')
title('Convergence history')
xlabel('Iteration number')
ylabel('||r||_2')



%CurlCurl0
% Warning: Matrix is close to singular or badly scaled. Results may be inaccurate. RCOND =  4.088716e-17. 
% > In blcLS (line 6)
%   In TestAllPreconditioners>@(x)defpreconditioner(Q2,R22,x)+blcLS(Pr,begin_in,end_in,x)
%   In iterapp (line 56)
%   In gmres (line 503)
%   In TestAllPreconditioners (line 154) 
% Warning: Matrix is close to singular or badly scaled. Results may be inaccurate. RCOND =  1.421495e-17.

