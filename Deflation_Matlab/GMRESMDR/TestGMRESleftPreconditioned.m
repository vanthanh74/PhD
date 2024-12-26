clear all
close all
clc

matrix='saylr4.mtx';  %sherman3 N=5005 nz=20033,   saylr4
A =  mmread(matrix);
fprintf('Size of A - %d\n', length(A));
fprintf('Non zeros of A - %d\n', nnz(A));

n = length(A);
npes =    8;  %number of subdomains nested disseciton npes=2^n
k=20;   % number of smallest singular values to deflate
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
V2_tilde = P_C*[-R11^-1*R12;eye(k)];  % size n x k
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

%scaling 
[A_kND_scaling,R00,C00] = equil_rar(A_kND);

% GMRES with left preconditioner

% Minv = (I-Q2Q2t +  Q2*R22^-1*Q2');
maxit = 2600;
A1 = A_kND_scaling;
H = zeros(maxit + 1, maxit);
V = zeros(n, maxit + 1);

x_0 = zeros(n, 1);
x_1 = x_0;
r_0 = preconditioner(Q2,R22,b-A1*x_0); %b-A1*x_0;%;% %preconditioner(Q2,R22,b-A1*x_0); %  % M^-1(b-Ax0)

% r_0 = rhs - lhs * x_0; % b-A*x_0

beta = norm(r_0);
V(:, 1) = r_0/beta;%(A1*rhs)/norm(A1*rhs);    % GMRES
% V(:, 1) = (A1*rhs)/norm(A1*rhs); %r_0/normr; %RRGMRES
% e_1 = zeros(maxit + 1, 1);
% e_1(1) = 1;
% c = beta * e_1;
resvec=[beta];

nb_subdomains = 8;
% % Partitioning using Metis
[p1, edgecut, s, begin_in, end_in] = DefineAlpha(A1, nb_subdomains);
% 
A2 = A1(p1, p1);
for i = 1 : nb_subdomains
    Pr{i} = A2(begin_in(i) : end_in(i), begin_in(i) : end_in(i));
end




% define M as a  function handle returning M\x
% aa = A1*V(:, j);
% bb = Q2'*aa;
% cc = R22\bb;
% dd = Q2*cc;
% ee = Q2*bb;
% w  = aa - ee + dd;
% M  = @preconditioner; 

for j = 1 : maxit       
%     w = A * blcLS(Pr, begin_in, end_in, V(:, j)); %blcLS  y(begin_in(i) : end_in(i), :) = Pr{i} \ x(begin_in(i) : end_in(i), :);    
%     w =  A_kND_scaling*blcLS(Pr, begin_in, end_in, V(:, j)); %blcLS  y(begin_in(i) : end_in(i), :) = Pr{i} \ x(begin_in(i) : end_in(i), :);    
    %   (I-Q2Q2t)*A_kND*x = (I-Q2Q2t)*b
% cond((I-Q2Q2t)*A_kND)
% cond((I-Q2Q2t + Q2*R22^-1*Q2')*A_kND)
%    1.6481e+18, 3.6977e+17

%     aa = A1*V(:, j);
%     bb = Q2'*aa;
%     cc = Q2*bb;
%     w  = aa - cc;
%-----------------------------------------------------------------------

%   (I-Q2Q2t + Q2*R22^-1*Q2')*A_kND*x = (I-Q2Q2t + Q2*R22^-1*Q2')*b
%     aa = w;
%     bb = Q2'*aa;
%     cc = R22\bb;
%     dd = Q2*cc;
%     ee = Q2*bb;
%     w  = aa - ee + dd;
%    w = preconditioner(Q2,R22,A2)*blcLS(Pr, begin_in, end_in, V(:, j)); 
%    w = A2*blcLS(Pr, begin_in, end_in, V(:, j)); 
     w = preconditioner(Q2,R22,A1*V(:, j)); 
%-----------------------------------------------------------------------
    
%     w =   A1* V(:, j);
    
    temp = V(:, 1 : j)' * w;
    w = w - V(:, 1 : j) * temp;
    temp1 = V(:, 1 : j)' * w;
    w = w - V(:, 1 : j) * temp1;                
    H(1 : j, j) = temp + temp1;
    H(j + 1, j) = norm(w);     
    V(:, j + 1) = w / H(j + 1, j);                 
    [QH, RH] = qr(H(1 : j + 1, 1 : j));       %    H = QR 
    e_1 = zeros(j + 1, 1);                    %    find y to minimize(||Hy - beta*e1||) , beta = ||b-Ax_0||
    e_1(1) = 1;                               %               ||QRy - beta*e1||
    c = beta * e_1;           % GMRES
%     c = [V(:,1 : j + 1 )]'*rhs ; % RRGMRES
    res = QH' * c;
    fprintf('Res(%d) = %e\n', j, abs(res(end)));
    resvec = [resvec, abs(res(end))];
%     y = RH(1 : j-1 , :) \ res(1 : j-1 );   
%     x_j = x_0 + V(:, 1 : j) * y;
%     x_0 = x_j;
%     r0 = preconditioner(Q2,R22,b-A1*x_0);
%     beta = norm(r0);
%     V(:,j+1) = r0/beta;
end
resvec=resvec';
%%
% k = maxdef;
% x_0 = zeros(n, 1);
y = RH \ res;       % solve y to minimize  ||Ry - Q'*beta*e1||
x_1 = x_0 + V(:, 1 : maxit) * y;
r_1 = V * (c - H * y);
            
% x_0 = zeros(length(b),1);
relres = norm(preconditioner(Q2,R22,b-A1*x_1))/norm(preconditioner(Q2,R22,b-A1*x_0))

% r_1 = V * (c - H * y);
% x_tilde = x_1;
% normX_Xexact = norm(x_exact-x_1)


%--------------------------------------------------------
%   (I-Q2Q2t + Q2*R22^-1*Q2')*A_kND*x = (I-Q2Q2t + Q2*R22^-1*Q2')*b
% M  = @(x) preconditioner(Q2,R22,x);% 
M1 = @(x) preconditioner(Q2,R22,x); 
% M2 = @(x) blcLS(Pr, begin_in, end_in, x);
% M  = @(x)  preconditioner(Q2,R22,A2)*blcLS(Pr, begin_in, end_in,x); 
%---------------------------------------------------------
% [x_gmres,flag,relresg,iter,resvec] = gmres(A1,b,[],1e-6,50); %
[x_gmres1,flag1,relres1,iter1,resvec1]=gmres(A1,b,[],1e-6,maxit,M1); % deflation preconditioner
relres1

norm(x_1-x_gmres1)
figure
plot(1:n,x_1,'b',1:n,x_gmres1,'.r')

