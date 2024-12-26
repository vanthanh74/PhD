clear all
close all
clc

matrix='sherman3.mtx';  %sherman3 N=5005 nz=20033,   saylr4
A =  mmread(matrix);
fprintf('Size of A - %d\n', length(A));
fprintf('Non zeros of A - %d\n', nnz(A));

n = length(A);
npes =    8;  %number of subdomains nested disseciton npes=2^n
k=20;   % number of smallest singular values to deflate
maxit = 25;
tol = 1e-8;

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
norm(A_ND*P_kND_mat-A_kND)
norm(A_ND-A_kND*P_kND_mat')

% A_kND*P_C = Q_final*R_final=  [Q1 Q2] *   [R11  R12]
%                                           [     R22]
% A_kND = full(A_kND);
[Q_final,R_final,P_C] = qr(full(A_kND),'matrix');   %   [Q,R,E] = QR(A) produces unitary Q, upper triangular R and a
%   permutation matrix E so that A*E = Q*R. The column permutation E is
%   chosen to reduce fill-in in R.
%----------------------------------------------------
% A_permutted = A_final;
normAPC_QR=normest(A_kND*P_C-Q_final*R_final)
Q1=Q_final(1:n,1:n-k);  % size n x n-k
Q2=Q_final(1:n,end-k+1:end); % size n x k
R11=R_final(1:n-k,1:n-k); % size n-k x n-k
R12=R_final(1:n-k,n-k+1:end); % size n-k x k
R22=R_final(end-k+1:end,end-k+1:end); % size k x k
normR22=norm(R22)
%%Approximate for the right null space of A
V2_tilde = P_C*[-R11^-1*R12;eye(k)];  % size n x k
Q2Q2t = Q2*Q2';
Q2tQ1 = Q2'*Q1;
PctPc=P_C'*P_C;
PcPct=P_C*P_C';
zerosIdentity=R22^-1*Q2'*A_kND*P_C;   % [0 I] size k x n
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
x_exact = A_kND\b;
% x_exact = gmres(A_kND,b,[],1e-6,800);
% sv qr deflated gmres

% to solve Ax = b, we solve A_kNDy=q
% x can be recovered by
% x =(I-V2_tilde*zerosIdentity*P_C')x_tilde + V2_tilde*R22^-1*Q2'*b
% where x_tilde is the approximate solution of (I-Q2Q2t)A_kNDx=(I-Q2Q2)b

%solve (I-Q2Q2t)A_kNDx=(I-Q2Q2)b

% lhs = (I-Q2Q2t)*A_kND;
% rhs = (I-Q2Q2t)*b;
% I (eye(size(Q2Q2t))


% %lhs
% bb = Q2'*A_kND;
% cc = Q2*bb;
% lhs=A_kND-cc;
% 
% %rhs
% dd = Q2'*b;
% ee = Q2*dd;
% rhs = b-ee;


% x_tilde = gmres(lhs,rhs,[],1e-6,400); % gmres(A,b,restart,tol,maxit,M1,M2,x,varargin)

% x_tilde = gmres((I-Q2Q2t+V2_tilde*(Q2'*A*V2_tilde)^-1*Q2')*A_kND,(I-Q2Q2t+V2_tilde*(Q2'*A*V2_tilde)^-1*Q2')*b);



%lhs
bb = Q2'*A_kND;
cc = R22\bb;
dd = Q2*cc;
ee = Q2*bb;
lhs0=A_kND - ee + dd;

%rhs
bb1 = Q2'*b;
cc1 = R22\bb1;
dd1 = Q2*cc1;
ee1 = Q2*bb1;
rhs0= b - ee1 + dd1;



% x_tilde = gmres((I-Q2Q2t   + Q2*R22^-1*Q2' )*A_kND,(I-Q2Q2t + Q2*R22^-1*Q2' )*b,[],1e-6,800); %gmres(A,b,restart,tol,maxit,M1,M2,x,varargin)
% x_tilde1 = gmres(lhs,rhs,[],1e-6,800); %gmres(A,b,restart,tol,maxit,M1,M2,x,varargin)
% % x_tilde3 = gmres(A_kND,b,[],1e-6,25,(I-Q2Q2t   + Q2*R22^-1*Q2' ));
% % x_tilde = ((I-Q2Q2t   + Q2*R22^-1*Q2' )*A_kND)\((I-Q2Q2t + Q2*R22^-1*Q2' )*b);
% norm_X_Xexact=norm(x_tilde-x_exact)
% norm_X_Xexact1=norm(x_tilde1-x_exact)
% norm_X_Xexact3=norm(x_tilde1-x_exact)
% x_tilde=pinv(full(A_kND))*b;

% x_tilde = cgs((I-Q2Q2t)*A_kND,(I-Q2Q2t)*b);

%--------------------------------------------------------------------
%BFGMRES
% x_tilde = BFGMRES((I-Q2Q2t    )*A_kND,(I-Q2Q2t)*b,tol,maxit);
% x_tilde_exact = ((I-Q2Q2t   + Q2*R22^-1*Q2' )*A_kND)\((I-Q2Q2t + Q2*R22^-1*Q2' )*b);
% norm(x_tilde-x_tilde_exact)
%--------------------------------------------------------------------
% %BFRRGMRES
% x_tilde = BFRRGMRES(lhs0,rhs0,tol,maxit);
% % x_tilde_exact = pinv(lhs)*rhs;
% x_tilde_gmres = gmres(lhs,rhs);
% norm(x_tilde-x_tilde_exact)
% norm(x_tilde-x_tilde_gmres)
%--------------------------------------------------------------------



%backslash matlab
% x_tilde= ((I-Q2Q2t)*A_kND)\((I-Q2Q2t)*b);
%--------------------------------------------------------------------

%pseudo inverse
% x_tilde= pinv((I-Q2Q2t)*A_kND)*((I-Q2Q2t)*b); 
%--------------------------------------------------------------------
% least squares method
% x_tilde= lsqr((I-Q2Q2t)*A_kND,(I-Q2Q2t)*b); 
%--------------------------------------------------------------------
% % GMRES

% % tol = 1e-6;
% % n = 500;
% % A = sprandn(n,n,0.1);
% % b = rand(n,1);
% % b=b/norm(b);
% % x_exact = A\b;
% % maxit = 25; 


% % (I-Q2Q2t)*A_kND*x = (I-Q2Q2t)*b
% % %lhs
% bb = Q2'*A_kND;
% cc = Q2*bb;
% lhs=A_kND-cc;
% 
% % %rhs
% dd = Q2'*b;
% ee = Q2*dd;
% rhs = b-ee;

% x1 = A_kND\b;
% x2 = ((I-Q2Q2t + Q2*R22^-1*Q2')*A_kND) \ ((I-Q2Q2t + Q2*R22^-1*Q2')*b);
% norm(x1-x2)

% (I-Q2Q2t + Q2*R22^-1*Q2')*A_kND*x = (I-Q2Q2t + Q2*R22^-1*Q2')*b
bb = Q2'*A_kND;
cc = R22\bb;
dd = Q2*cc;
ee = Q2*bb;
lhs=A_kND - ee + dd;
% cond(lhs)
% cond(((I-Q2Q2t + Q2*R22^-1*Q2')*A_kND))
%rhs
bb1 = Q2'*b;
cc1 = R22\bb1;
dd1 = Q2*cc1;
ee1 = Q2*bb1;
rhs= b - ee1 + dd1;

% x3 = lhs\rhs;
% norm(x1-x3)
% norm(x2-x3)
maxit = 400;
A1 = A_kND;
H = zeros(maxit + 1, maxit);
V = zeros(n, maxit + 1);

x_0 = rand(n, 1);
x_1 = x_0;
    
r_0 = b - A1 * x_0; % b-A*x_0

normr = norm(r_0);
V(:, 1) = r_0/normr;
e_1 = zeros(maxit + 1, 1);
e_1(1) = 1;
c = normr * e_1;
resvec=[normr];

% nb_subdomains = 8;
% % Partitioning using Metis
% [p1, edgecut, s, begin_in, end_in] = DefineAlpha(A1, nb_subdomains);
% 
% A2 = A1(p1, p1);
% for i = 1 : nb_subdomains
%     Pr{i} = A2(begin_in(i) : end_in(i), begin_in(i) : end_in(i));
% end
for j = 1 : maxit       
%     w = A * blcLS(Pr, begin_in, end_in, V(:, j)); %blcLS  y(begin_in(i) : end_in(i), :) = Pr{i} \ x(begin_in(i) : end_in(i), :);    
%   (I-Q2Q2t + Q2*R22^-1*Q2')*A_kND*x = (I-Q2Q2t + Q2*R22^-1*Q2')*b
    aa = A1*V(:, j);
    bb = Q2'*aa;
    cc = R22\bb;
    dd = Q2*cc;
    ee = Q2*bb;
    w  = aa - ee + dd;


%     w =   A1* V(:, j);
    
    temp = V(:, 1 : j)' * w;
    w = w - V(:, 1 : j) * temp;
    temp1 = V(:, 1 : j)' * w;
    w = w - V(:, 1 : j) * temp1;                
    H(1 : j, j) = temp + temp1;
    H(j + 1, j) = norm(w);     
    V(:, j + 1) = w / H(j + 1, j);                 
    [QH, RH] = qr(H(1 : j + 1, 1 : j));   %    H = QR 
    e_1 = zeros(j + 1, 1);                 % find y to minimize(||Hy - beta*e1||)
    e_1(1) = 1;                               %               ||QRy - beta*e1||
    c = normr * e_1;                             
    res = QH' * c;
    fprintf('Res(%d) = %e\n', j, abs(res(end)));
    resvec = [resvec, abs(res(end))];
end

% k = maxdef;
y = RH(1 : end - 1, :) \ res(1 : end - 1);       % solve y to minimize  ||Ry - Q'beta*e1||
x_1 = x_0 + V(:, 1 : maxit) * y;
% r_1 = V * (c - H * y);
% x_tilde = x_1;
normX_Xexact = norm(x_exact-x_1)
[x_gmres,flag,relres,iter,resvec1]=gmres(A_kND,b,[],1e-6,maxit,[],[],x_0);
normXgmres_X = norm(x_1-x_gmres)

figure
plot(1:length(resvec1),resvec1,'b')
hold on
plot(1:length(resvec),resvec,'.r')
%-------------------------------------------------------

%recover x
% x = (I-V2_tilde*zerosIdentity*P_C')*x_tilde  +  V2_tilde*pinv(R22)*Q2'*b;
% x = x_tilde;
% aa=(I-Q2Q2t)*A_kND;
% det(aa);
x = x_1;
x_gmres = gmres(A_kND,b);
norm_X_Xexact=norm(x-x_exact)
norm_X_Xgmres=norm(x-x_gmres)
norm_Xgmres_Xexact=norm(x_gmres-x_exact)
%-------------------------------------------------------
figure
plot(x)
legend('x')

figure
plot(x_exact)
legend('x_{exact}')

figure
plot(x)
hold on
plot(x_exact)
legend('x','x_{exact}')





% % svd deflated gmres
% maxit=100
% A = rand(50,50);
% b=rand(50,1);
% b=b/norm(b);
% n=length(A)
% H = zeros(maxit + 1, maxit);
% V = zeros(n, maxit + 1);
% 
% x_0 = zeros(n, 1);
% x_1 = x_0;
% r_0 = b - A * x_0;
% 
% normr = norm(r_0);
% V(:, 1) = r_0/normr;
% e_1 = zeros(maxit + 1, 1);
% e_1(1) = 1;
% c = normr * e_1;
% resvec=[];
% for j = 1 : maxit
% %     w = A * blcLS(Pr, begin_in, end_in, V(:, j)); %blcLS  y(begin_in(i) : end_in(i), :) = Pr{i} \ x(begin_in(i) : end_in(i), :);
%     w = A * V(:, j);
%     temp = V(:, 1 : j)' * w;
%     w = w - V(:, 1 : j) * temp;
%     temp1 = V(:, 1 : j)' * w;
%     w = w - V(:, 1 : j) * temp1;
%     H(1 : j, j) = temp + temp1;
%     H(j + 1, j) = norm(w);
%     V(:, j + 1) = w / H(j + 1, j);
%     [QH, RH] = qr(H(1 : j + 1, 1 : j));
%     e_1 = zeros(j + 1, 1);
%     e_1(1) = 1;
%     c = normr * e_1;
%     res = QH' * c;
%     fprintf('Res(%d) = %e\n', j, abs(res(end)));
%     resvec = [resvec, abs(res(end))];
% end
% 
% % k = maxdef;
% y = RH(1 : end - 1, :) \ res(1 : end - 1);
% x_1 = x_0 + V(:, 1 : maxit) * y;
% r_1 = V * (c - H * y);
% 
% norm(A\b-x_1)
% 
% %----------------------------------
% % svd deflation
% [AP_k, Sigma_k, P_k] = svd(H,0);
% P_k = P_k(:, end - k + 1 : end);
% Y_k = V(:, 1 : maxit) * P_k;
% % [Q, R] = qr(A * blcLS(Pr, begin_in, end_in, Y_k), 0);
% [Q, R] = qr(A * Y_k, 0);
% C_k = Q;
% U_k = Y_k / R;
% %----------------------------------
% 
% figure
%  plot(log10(resvec));
%  figure
%  plot(x_1)
