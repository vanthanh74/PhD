% function x_k = BFGMRES(A,b,tol,maxit)

clear all
close all
clc
%%
% matrix='sherman3.mtx';  %sherman3 N=5005 nz=20033,   saylr4
% A =  mmread(matrix);
% fprintf('Size of A - %d\n', length(A));
% fprintf('Non zeros of A - %d\n', nnz(A));
% 
% n = length(A);
% npes =    8;  %number of subdomains nested disseciton npes=2^n
% k=20;   % number of smallest singular values to deflate
% maxit = 25;
% tol = 1e-1;
% N = 0.5*tol;
% 
% %A_ND = A(p,p),     [p,ip,sizes] = metismex('NodeNDP',A,npes); %
% % matrix A_kND obtain from QRCP with nested dissection with the last k
% % columns corresponding to the k smallest singular values
% [A_ND,A_kND,P_kND]= NDQRCPSingularValuesApproximationPk(A,npes,k);



%%
% % breakdown-free GMRES
% N = 0.5;
% tol = 1e-8;
% maxit = 25;
% n = 5000;
% A = sprandn(n,n,0.0001);
% b = rand(n,1);
% b=b/norm(b);

%  test = 'test.mat';
% A = (I-Q2Q2t)*A_kND;
% b = ((I-Q2Q2t)*b);
% save(test,'A','b')
load('test.mat')
%%
% A=A_kND;
maxit = 25;
n = length(A);
tol = 1e-1;
N = 0.5*tol;
x_exact = A\b;
% x_exact = pinv(full(A))*b;
 %%
x_0 = zeros(n,1);
x_k = x_0;
% % x_0=x_0/norm(x_0);
% r_0 = b - A * x_0;
% normr = norm(r_0);

v1 = (A*b)/norm(A*b);
p = 0;
V = [v1];
U = [];
H_hat = [];
G = [];
resvec=[];
    i=0;
r_0=b;
         
for k=1:maxit
    %     k=k+1
    w = A*V(:,k);
    hk = V(:, 1 : k)' * w;
    if size(U,1)~=0
        gk = U(:, 1 : p)' * w;
        w = w - V(:, 1 : k)*hk - U(:, 1 : p)*gk;
    else
        w = w - V(:, 1 : k)*hk ;
    end
    
    hk1 = norm(w);
    H_hat_temp = H_hat;
    H_hat = [H_hat_temp hk; zeros(1,k-1) hk1];
    
    G_temp = G;

    if cond(H_hat) > 10^(2*p)/tol %|| norm( x_k-x_0) < N*norm(x_k)
        U(:,p+1) = V(:,k);

        G_temp = [G_temp; H_hat_temp(k,:)];
        H_hat_temp(k,:) = 0;
        
%         233 190 136 126 126 156 
        while cond(H_hat) > 10^(2*p)/tol 
                    i=i+1
            v_hat = rand(length(b),1);%A'*r_k;%
            temp1 = V(:, 1 : k-1)' * v_hat;
            temp2 = U(:,1 : p + 1)' * v_hat;
            v_hat = v_hat - V(:, 1 : k-1)*temp1- U(:,1 : p + 1)*temp2;
            v_hat = v_hat/norm(v_hat);                       
            V(:,k) = v_hat;
                
            w = A*V(:,k);
            hk = V(:, 1 : k)' * w;
            gk = U(:, 1 : p + 1)' * w;
            w = w - V(:, 1 : k)*hk - U(:, 1 : p + 1)*gk;
            
            hk1 = norm(w);
            H_hat = [H_hat_temp hk; zeros(1,k-1) hk1];
%             cond(H_hat(1:k,1:k) )
            cond(H_hat)
        
            %         if cond(H_hat(:,1:k)) >  10^(2*p)/tol
            %             v_hat = rand(length(b),1);
            %             temp1 = V(:, 1 : k-1)' * v_hat;
            %             temp2 = U(:,1 : p + 1)' * v_hat;
            %             v_hat = v_hat - V(:, 1 : k-1)*temp1- U(:,1 : p + 1)*temp2;
            %             v_hat = v_hat/norm(v_hat);
            %             V(:,k) = v_hat;
            %         end
            %             cond(H_hat(:,1:k))
            %              p = p+1;
            %              U(:,p+1) = V(:,k);
        end
        p = p+1                 
    end
    
    vk1 = w/hk1;
    V(:,  k+1) = vk1;
    
    if p > 0
        G = [G_temp,gk];
    else
        G = G_temp;
    end
    
    H_hatG = [H_hat;G];
    [QH, RH] = qr(H_hatG(1 : k + 1 + p, 1 : k));
%     e_1 = zeros(k + 1, 1);
%     e_1(1) = 1;
    c = [V(:,1 : k + 1 ),U(:,1:p)]'*b ;
    res = QH' * c;
    fprintf('Res(%d) = %e\n', k, abs(res(end)));
    resvec = [resvec, abs(res(end))];
    
    y = RH(1 : end - 1, :) \ res(1 : end - 1);
    x_0 = x_k;
    x_k =  V(:, 1 : k) * y;
    r_k = b-A*x_k;
    r0=r_k;
%     x_k = x_0 +  V(:, 1 : k) * y;

    %     r_1 = V * (c - H_hatG * y);
    %     normr = norm(r_1);
    
end
% y = RH(1 : end - 1, :) \ res(1 : end - 1);
% x_k =  V(:, 1 : maxit) * y;
% r_1 = V * (c - H_hatG * y);

x_gmres = gmres(A,b);
normX_gmres_X = norm(x_gmres-x_k)/norm(x_k)
normX_Xexact = norm(x_exact-x_k)/norm(x_k)
normXexact_X_gmres = norm(x_gmres-x_exact)/norm(x_gmres)
% hold on
% plot(x_k)
% plot(x_gmres)
% V(:,2)'*V(:,8)
% V(:,1)'*v_hat
