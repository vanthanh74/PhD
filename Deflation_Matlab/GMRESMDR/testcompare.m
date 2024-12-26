clear all
close all
clc

% matrix='saylr4.mtx';  %sherman3 N=5005 nz=20033,  bcsstm25 bcsstk35  saylr4 lshp3466 CurlCurl_0 crystk03.mtx ct20stif.mtx bcsstm39.mtx

matrix='saylr4.mtx';  %sherman3 N=5005 nz=20033,   saylr4, lshp3466
A =  mmread(matrix);
% e20r0100 e05r0400 e30r0100
% figure
% spy(A)
% fprintf('Size of A - %d\n', length(A));
% fprintf('Non zeros of A - %d\n', nnz(A));

n = length(A);
% n = 500; % size of the matrix A
npes =    8;  %number of subdomains nested disseciton npes=2^n
k_df= 50   % number of smallest singular values to deflatef % rank of R22
k = n-k_df; % rank of A11


%%-------------------------------------------------------------------------
% % generate matrix with desired singular values
% % n = 500; % size of the matrix A
% S = rand(n,1); % singular values of A=USV^T
% S = sort(S,'descend');
% S = diag(S); % singular values of A=USV^T
% S(end-k_df+1:end) = 1e-12*S(end-k_df+1:end);
% 
% figure
% plot(1:n,S,'.')
% 
% % U = rand(n);
% U = sprand(n,n,0.001);
% [Uu,Ru] = qr(U);
% U=Uu;
% % V = rand(n);
% V = sprand(n,n,0.001);
% [Vv,Rv] = qr(U);
% V=Vv;
% 
% %recover A
% A = U*S*V';
% A=sparse(A);
% s=svd(full(A));
% norm(S-s)
% figure
% spy(A)
% % s=svd(A);
% % Ss=sort(diag(S),'descend')
% % % -----------------------------------------------------------------------
% A=sparse(A);
fprintf('Size of A - %d\n', length(A));
fprintf('Non zeros of A - %d\n', nnz(A));



%scaling 
[A,R00,C00] = equil_rar(A);


% n = 500;9
% A = sprandn(n,n,0.01);
% figure
% spy(A)


% k = 20;
% maxit = 25;
tol = 1e-6;

%A_ND = A(p,p),     [p,ip,sizes] = metismex('NodeNDP',A,npes); %
% matrix A_kND obtain from QRCP with nested dissection with the last k
% columns corresponding to the k smallest singular values
[A_ND,A_kND,P_kND,p,last_k_columns_matrix,last_k_columns_matrix_temp]= NDQRCPSingularValuesApproximationPk(A,npes,k,k_df);



% A_kND = A_ND*P_kND
% A_ND = A(p,p) = P'*A*P
% figure
% spy(A_ND)
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



I=eye(n,n);
P_ND = I(:,p); %P_ND
Pi = I(:,P_kND); % Pi
% normest(A_ND*P_kND_mat-A_kND)
% norm(A_ND-A_kND*P_kND_mat')
normest(A_kND-P_ND'*A*P_ND*Pi)
% A_kND*P_C = Q_final*R_final=  [Q1 Q2] *   [R11  R12]
%                                           [     R22]


%----------------------------------------------------
% % k-way partition
% nb_subdomains = 8;
% % % Partitioning using Metis
% 
% [p1, edgecut, s, begin_in, end_in] = DefineAlpha(A, nb_subdomains);
% P_kw = I(:,p1);    %  P_kw
% A = A(p1, p1); % A = P_kw'*A*P_kw
% for i = 1 : nb_subdomains
%     Pr{i} = A(begin_in(i) : end_in(i), begin_in(i) : end_in(i));      
% end
% 
% P=A0(begin_in : end_in, begin_in : end_in);    
%----------------------------------------------------
% figure
% spy(A(p,p))

% figure spy( Pr{2})
% A_kND = full(A_kND);
% [Q_final,R_final] = qr(full(A_kND),'matrix'); 
% do strong RRQR for A_kND
% k = 40;

% f = 2;
% [Q_final, R_final, P_final] = sRRQR(full(A_kND), f, 'rank', k);

[Q_final, R_final] = qr(full(A_kND),'matrix');


% [Q_final,R_final] = qr(full(A_kND),'matrix');      
% P_kND_mat=P_C;
%----------------------------------------------------
% A_permutted = A_final;
% normAPC_QR=normest(P_kND*P_C-Q_final*R_final)
% Q1=Q_final(1:n,1:n-k);  % size n x n-k
Q2=Q_final(1:n,end-k_df+1:end); % size n x k
R11=R_final(1:n-k_df,1:n-k_df); % size n-k x n-k
R12=R_final(1:n-k_df,n-k_df+1:end); % size n-k x k
R22=R_final(end-k_df+1:end,end-k_df+1:end); % size k x k
normR22=norm(R22)
%%Approximate for the right null space of A
V2_tilde = P_kND_mat*[-R11^-1*R12;eye(k_df)];  % size n x k
% Q2Q2t = Q2*Q2';
% Q2tQ1 = Q2'*Q1;
% PctPc=P_C'*P_C;
% PcPct=P_C*P_C';
zerosIdentity=R22\(Q2'*A_ND*P_kND_mat);   % [0 I] size k x n
% aaa=V2_tilde*V2_tilde^-1;

%verify
% Q22 = A_kND*V2_tilde*R22^-1;
% V2tilde=A_kND^-1*(Q2*R22);
% AV2tilde=A_ND*V2_tilde;
% Q2R22=Q2*R22;
% norm(AV2tilde-Q2R22)
% R22_verify = Q2'*A_ND*V2_tilde;
% norm(R22-R22_verify)

%------------------------------------------------------------------
% M^-1 * A_kND
% MinvA = (I-Q2*Q2'+V2_tilde*(R22)^-1*Q2')*A_kND;
% rhsMinvA=[Q1 V2_tilde]*[R11 R12; zeros(k,n-k) eye(k)];
% normMinvA=norm(MinvA-rhsMinvA)
% sA=svd(full(A_kND));
% % sA_kND = sA;
% % sA1=svd(full(A));
% % sA2=svd(full(A_ND));
% % sk1=sA1(end-k:end);
% % sk2=sA2(end-k:end);
% 
% s=svd(MinvA);
% sA_k = sA(end-k:end);
% s_k = s(end-k:end);
% 
% 
% % aa=R22^-1*Q2'*A_kND - zerosIdentity*P_C';
% % aa=V2_tilde*zerosIdentity*P_C';
% % bb=aa*P_C';
% %------------------------------------------------------------------
% % M^-1 * A_kND
% MinvA1 = (I-Q2*Q2' +  Q2*R22^-1*Q2')*A_kND;
% rhsMinvA1=[Q1 Q2]*[R11 R12; zeros(k,n-k) eye(k)];
% normMinvA1=norm(MinvA1-rhsMinvA1)
% s1=svd(MinvA1);
% s_k1 = s1(end-k:end);

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

% x_exact = gmres(A_kND,b,20,1e-1)
% x_exact = pinv(full(A_kND))*b;
x_exact = A_ND\b;


%x_gmres
[x_gmres,flag_gmres,relres_gmres,iter_gmres,resvec_gmres] = gmres(A_ND,b,[],1e-6,length(b));  %



% %-------------------------------------------------------
%svd deflation
[U,Sigma,V] = svd(full(A_ND));
U2 = U(:,end-k_df+1:end);
Sigma2 = Sigma(end-k_df+1:end,end-k_df+1:end);
V2 = V(:,end-k_df+1:end);  

%solve (I-U2U2')Ax = (I-U2U2')b for x_tilde
% I_k = I(1:k,1:k);
% x_tilde_svd = ((I-U2*U2')*A_ND)\((I-U2*U2')*b);
[x_tilde_svd,flag_svd,relres_svd,iter_svd,resvec_svd] = gmres((I-U2*U2')*A_ND,(I-U2*U2')*b,[],1e-6,length(b));  %
x_svd = (I-V2*V2')*x_tilde_svd + V2*(Sigma2\U2')*b;
norm(x_gmres-x_svd)   % sherman3 ans = 5.3148e+09, saylr4 5.4410e-09, lshp3466 1.8390e-11, msc01440 1.6074e-11

% %-------------------------------------------------------
%qr deflation
%solve (I-Q2Q2')Ax = (I-Q2Q2')b for x_tilde
% x_tilde_qr = ((I-Q2*Q2')*A_ND)\((I-Q2*Q2')*b);


[x_tilde_qr,flag_qr,relres_qr,iter_qr,resvec_qr] = gmres((I-Q2*Q2')*A_ND,(I-Q2*Q2')*b,[],1e-6,length(b));  %
x_qr = (I-V2_tilde*zerosIdentity*P_kND_mat')*x_tilde_qr  +  V2_tilde*(R22\Q2')*b;
norm(x_gmres-x_qr)     % sherman3 ans = 0.0112,  saylr4 7.5101e-08, lshp3466  1.1586e-11, msc01440 2.0352e-12


figure
plot(x_exact)
hold on
plot(x_qr)

figure
semilogy(resvec_qr)
hold on
semilogy(resvec_svd)
semilogy(resvec_gmres)
legend('Strong RRQR deflated GMRES','SVD deflated GMRES','GMRES')
title('Strong RRQR and SVD deflation')
xlabel('Iteration number')
ylabel('||r||_2')
% xlim([0 5262])
xlim([0 1000])
