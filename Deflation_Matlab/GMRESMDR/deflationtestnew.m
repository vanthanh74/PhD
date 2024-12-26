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
[A0,R00,C00] = equil_rar(A);
% svdAscaling=svd(full(A));
% A0=A;

n = length(A);
npes =    8;  %number of subdomains nested disseciton npes=2^n
k_df= 50   % number of smallest singular values to deflatef % rank of R22
k = n-k_df; % rank of A11
% k = 20;
% maxit = 25;
tol = 1e-6;

%A_ND = A(p,p),     [p,ip,sizes] = metismex('NodeNDP',A,npes); %
% matrix A_kND obtain from QRCP with nested dissection with the last k
% columns corresponding to the k smallest singular values
[A_ND,A_kND,P_kND,p,last_k_columns_matrix,last_k_columns_matrix_temp]= NDQRCPSingularValuesApproximationPk(A0,npes,k,k_df);


%----------------------------------------------------
%verify
% A = randi([1 8],8)
% I=eye(size(A));
% p=[1 3 5 7 2 4 6 8];
% I(:,p)
% A*I(:,p)

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
%     
%     
% end
%----------------------------------------------------
% figure
% spy(A(p,p))

% A_kND = full(A_kND);
[Q_final,R_final] = qr(full(A_kND),'matrix'); 
% do strong RRQR for A_kND
% k = 40;

f = 2;
% [Q_final, R_final, P_final] = sRRQR(full(A_kND), f, 'rank', k);


% k_df=k;
% normest(A_kND-Q_final*R_final)

% Q_tilde = P_kw'*P_ND*Q_final;
Q_tilde = Q_final;


% R_tilde = R_final*Pi'*P_ND'*P_kw;
% Q_tilde1 = P_kw'*P_ND'*Q_final;

% A_permutted = A_final;
P_C = Pi;
% normAPC_QR=normest(A_kND*P_C-Q_final*R_final)
% Q1=Q_tilde(1:n,1:n-k);  % size n x n-k
Q2=Q_tilde(1:n,end-k_df+1:end); % size n x k
% R11=R_final(1:n-k,1:n-k); % size n-k x n-k
% R12=R_final(1:n-k,n-k+1:end); % size n-k x k
R22=R_final(end-k_df+1:end,end-k_df+1:end); % size k x k
normR22=normest(R22)




%----------------------------------------------------


% figure
% subplot(121)
% spy(R_final)
% subplot(122)
% spy(R22)
% figure
% spy(R22)
% 
% 
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

% %svd(A)
% sA=svd(full(A));
% % [u,sA,v]=svd(full(A));
% k_smallest_sv_svd = sA(end-k+1:end);
% 
% % sA_kND = svd(full(A_kND));
% % sA_kND_k =sA_kND(end-k+1:end);
% 
% %svd(R22)
% sR=svd(R_final);
% sR_k = sR(end-k+1:end);
%  
% sR22=svd(R22);
% sR22_k = sR22(end-k+1:end);
% 
% [k_smallest_sv_svd  sR_k sR22_k]
% 
% figure
% spy(R_22)
% 
% figure
% semilogy(1:k,k_smallest_sv_svd,'-db',1:k,sR22_k,'-^k')
% % legend('SV by SVD(A)','Approximate SV by QRCP on A')
% legend_string={'SV by SVD(A)','Approximate SV by QRCP on A'};


% figure
% spy(R_final)
% figure
% spy(R_tilde)
% s=svd(full(R_final));
% s1=svd(full(R_tilde));
% normest(s-s1)

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
sA=svd(full(A0));
% sA=svd(full(A_kND));
% % [u,sA,v]=svd(full(A));
sA_k = sA(end-k_df+1:end);
% sA_k = svds(A,k,'smallest') ;
   
% last_6_column_A = full(A(:,end-6+1:end))
% last_6_column_selected_A = full(A(:,P_kND(end-6+1:end)))
% norm(last_6_column_A-last_6_column_selected_A)
% sA_k0 = diag(sA(end-k+1:end,end-k+1:end));
% sA_kND=svd(full(A_kND));
% sA_kkND = sA_kND(end-k+1:end);
%------------------------------------------------------------------

%svd(R_final)
sR=svd(full(R_final));
sR_k = sR(end-k_df+1:end);
% sR_k = svds(R_final,k,'smallest') ;
%------------------------------------------------------------------

% % R_22
sR22=svd(full(R22));
sR22_k = sR22(end-k_df+1:end);
% sR22_k = svds(R22,k,'smallest') ;
%------------------------------------------------------------------
% % M^-1 * A
% matrix='CurlCurl_0.mtx'; 
% A0 =  mmread(matrix);
MinvA = defpreconditioner1(Q2,R22,P_ND,A0);%;%(I-Q2*Q2' +  Q2*(R22^-1)*Q2')*A;%

% %A-DEF2 preconditioner
% QQ = Q2*(Q2'*A*Q2)^-1*Q2'; %Qfunc
% P = I - A*QQ;
% MinvA = P'*defpreconditioner(Q2,R22,A) + QQ;%;%(I-Q2*Q2' +  Q2*(R22^-1)*Q2')*A;%


% MinvA00 =  defpreconditioner(Q2,R22,A);
% bb = Q2'*A;
% cc = R22^-1*bb;
% dd = Q2*cc;
% ee = Q2*bb;
% MinvA0  = A - ee + dd;
% normest(MinvA00-MinvA0)
% normest(MinvA-MinvA0)
% normest(MinvA00-MinvA)

s=svd(full(MinvA));
s_k = s(end-k_df+1:end);
% s_k = svds(MinvA,k,'smallest') ;

% s0=svd(MinvA0);
% s_k0 = s0(end-k+1:end);

%------------------------------------------------------------------
% % (M_1^-1 + M_2^-1) * A
% aa=rand(length(A),1);
% blJacobi=blcLS(Pr, begin_in, end_in,aa);
% diagblJacobi = diag(blJacobi);
% n = length(begin_in);
% for i = 1 : n
%          y(begin_in(i) : end_in(i), :) = full(Pr{i})* A(begin_in(i) : end_in(i), :);%Pr{i}\ x(begin_in(i) : end_in(i), :);%Pr{i}\ x(begin_in(i) : end_in(i), :);
% end




% (I-Q2*Q2' +  Q2*(R22^-1)*Q2')*A
% MinvA1 = MinvA + blcLS(Pr, begin_in, end_in,A);%defpreconditioner(Q2,R22,A) + blcLS(Pr, begin_in, end_in,A);
% MinvA1 = defpreconditioner(Q2,R22,blcLS(Pr, begin_in, end_in, A)); % (M1^-1M2^-1)A
% MinvA1 = defpreconditioner1(Q2,R22,P_ND,A0) + blcLS(Pr, begin_in, end_in,A) - blcLS(Pr, begin_in, end_in,A)*defpreconditioner1(Q2,R22,P_ND,A);


% s1=svd(full(MinvA1));
% s_k1 = s1(end-k_df+1:end);
% s_k1 = svds(MinvA1,k,'smallest') ;
%------------------------------------------------------------------
% %  M_2^-1 * A
% MinvA2 =  blcLS(Pr, begin_in, end_in,A);
% s2=svd(full(MinvA2));
% s_k2 = s2(end-k_df+1:end);
%  s_k2 = svds(MinvA,k,'smallest') ;  
dd=1
compare_mat=[sA_k sR_k sR22_k s_k  ]
% figure
% subplot(131)
% spy(MinvA) % def
% subplot(132)
% spy(MinvA1) % mixed
% subplot(133)
% spy(MinvA2) % block Jacobi
