clear all; clc; close all;
format short e;


% adapt to span the deflation subspace Y_k = V_k*P_k
% save results
write_results = 0;

% addpath('matrix');
% The test matrices
% Pr_BIGP0='matrix/BIGP1_METRIC1_ELX_TSTEP_114_NEWTONIT_0_Pressure.mtx';
% Pr_BIGP1='matrix/BIGP1_METRIC1_ELX_TSTEP_114_NEWTONIT_1_Pressure.mtx';

% Pr_BIGP0='fidap001.mtx';%N=216
% Pr_BIGP0='fidap003.mtx';%N=1821
% Pr_BIGP0='fidap021.mtx';%N=656
% Pr_BIGP0='fidap024.mtx';%N=2283 nz=47897
% Pr_BIGP0='g7jac040sc.mtx';
% Pr_BIGP0='cavity01.mtx';% N=317
% Pr_BIGP0='dw4096.mtx';%N=8192 nz=41746
% Pr_BIGP0='msc01440.mtx'; %N=1440
% Pr_BIGP0='msc23052.mtx';%N=23052
% Pr_BIGP0='saylr4.mtx';  % N=3564 nz=22316
Pr_BIGP0='sherman3.mtx';  % N=5005 nz=20033
% Pr_BIGP0='sherman5.mtx';  % N=3312 nz=20793
% Pr_BIGP0='orsreg_1.mtx';  % N=2205 nz=14133



% Pr_BIGP1='fidap021.mtx';
% matrices = {Pr_BIGCO0, Pr_BIGCO0, Pr_BIGCO0, Pr_BIGCO1, Pr_BIGCO1, Pr_BIGCO1, Pr_BIGCO2, Pr_BIGCO2, Pr_BIGCO2};
% matrices = {Pr_BIGP0, Pr_BIGP0, Pr_BIGP0, Pr_BIGP0, Pr_BIGP0, ...
%             Pr_BIGP1, Pr_BIGP1, Pr_BIGP1, Pr_BIGP1, Pr_BIGP1, ...
%             Pr_BIGP2, Pr_BIGP2, Pr_BIGP2, Pr_BIGP2, Pr_BIGP2};

% matrices = {Pr_BIGP0, Pr_BIGP0, Pr_BIGP1, Pr_BIGP1};
matrices = {Pr_BIGP0}
% matrices = {Pr_BIGCO0, Pr_BIGCO0, Pr_BIGCO0, Pr_BIGCO0, Pr_BIGCO0...
%             Pr_BIGCO0, Pr_BIGCO0, Pr_BIGCO0, Pr_BIGCO0, Pr_BIGCO0...
%             Pr_BIGCO1, Pr_BIGCO1, Pr_BIGCO1, Pr_BIGCO1, Pr_BIGCO1...
%             Pr_BIGCO1, Pr_BIGCO1, Pr_BIGCO1, Pr_BIGCO1, Pr_BIGCO1...
%             Pr_BIGCO2, Pr_BIGCO2, Pr_BIGCO2, Pr_BIGCO2, Pr_BIGCO2...
%             Pr_BIGCO2, Pr_BIGCO2, Pr_BIGCO2, Pr_BIGCO2, Pr_BIGCO2};



% matrices = {Pr_BIGP0, Pr_BIGP0, Pr_BIGP0, Pr_BIGP0, Pr_BIGP0, ...
%             Pr_BIGP1, Pr_BIGP1, Pr_BIGP1, Pr_BIGP1, Pr_BIGP1, ...
%             Pr_BIGP2, Pr_BIGP2, Pr_BIGP2, Pr_BIGP2, Pr_BIGP2};
% matrices = {Pr_BIGP0, Pr_BIGP0, Pr_BIGP0, Pr_BIGP0, Pr_BIGP0, Pr_BIGP0, Pr_BIGP0, Pr_BIGP0, Pr_BIGP0, Pr_BIGP0, Pr_BIGP0};
% matrices = {BIGP52, BIGP52, BIGP52};
%   matrices = {sky2d, sky2d, sky2d};
%   matrices = {sky3d, sky3d, sky3d};
%  matrices = {ela25, ela25, ela25};
% matrices = {ela12, ela12};
% matrices = {LS3, LS3, LS3};
% matrices = {SPE10};


% Load Matrix
% kappa = 1;
% N = 100;
% [A,b,x] = gravity(n);
% [A,b,x] = heat(N,kappa);

matrix = matrices{1};
% matrix = A;
% A1 = A;
A1 = mmread(matrix);
%  [A1,R,C] = equil_rar(A1);

fprintf('Size of A - %d\n', length(A1));
fprintf('Non zeros of A - %d\n', nnz(A1));
numeq = length(A1);
npes =    8;  %number of subdomains nested disseciton npes=2^n
nb_subdomains = 8;%2*npes-1;   % number of blocks
tol = 1e-1;
tol_1 = 1e-1;
tol_2 = 1e-1;
iter_max = 2500;
%%%

rhs = rand(length(A1), 1);
rhs = rhs/norm(rhs);

% Safecheck because if nb_subdomains == 1 then Matlab crashes...
if nb_subdomains == 1
    fprintf('ERROR: the number of subdomains should be > 1!!!\n');
    return
end

% Partitioning using Metis
% [p1, begin_in, end_in, sizes] = DefineAlphaND(A1, npes);
% [p1, edgecut, s, begin_in, end_in] = DefineAlpha(A1, nb_subdomains);

% 1.    p1 the permutation vector of metismex(kway)
% 2.    edgecut the number of edges that are split between different
%       domains(output of metismex)
% 3.    s contains the sizes of each of the subdomains
% 4.    beginIn the start indices of the subdomains
% 5.    endIn the end indices of the subdomains

% [p,iperm] = metismex('EdgeND',A1);

%[perm,iperm] = metismex('EdgeND',A,options)
%[perm,iperm,sizes] = metismex('NodeNDP',A,npes,options)

% A = A1(p,p);
% A = A1(p1,p1);
n = length(A1);
%
%
% figure
% subplot(121)
% spy(A1(p1,p1))
% subplot(122)
% spy(A1(p1,p1))
%
% figure
% spy(A)





% figure
% spy(A(begin_in(1):end_in(1),begin_in(1):end_in(1)))




% Block Jacobi
Pr = cell(nb_subdomains, 1);

% Parameters of GMRES
maxit = 25;%50;%
maxdef = 5;%10;%
maxcycle = 5;
k = 0;
rng(1);
% Matrices of Krylov
H = zeros(maxit + 1, maxit);
V = zeros(n, maxit + 1);
resvec = cell(length(matrices), 3);
resvecND = cell(length(matrices), 3);

% Method of recycling
adaptive_deflation = 0;
reduction = 0;
methodDef = {'svd','reig', 'heig'};
methodDef = {'svd'};
for methodInd = 1 : length(methodDef)
    method = methodDef{methodInd};
    
    for matInd = 1 : length(matrices)
        cycle = 1;
        %         if (matInd == 1)
        %             tol = tol_1;
        %             if( adaptive_deflation == 1 )
        %                method = 'reig';
        %             end
        %         else
        %             tol = tol_2;
        %             method = methodDef{methodInd};
        %         end
        resvec{matInd, methodInd} = [];
        resvecND{matInd, methodInd} = [];
        A1 = mmread(matrices{matInd});
        [A1,R,C] = equil_rar(A1);
        A=A1;
        %         A = A1(p1, p1);
        %         figure
        %         spy(A)
        %         for i = 1 : length(begin_in)
        %             if begin_in(i)*end_in(i)~=0
        %                 Pr{i} = A(begin_in(i) : end_in(i), begin_in(i) : end_in(i));
        %             end
        %         end
        %         for ii=1:length(Pr)
        %             if size(Pr(ii),1)*size(Pr(ii),2) == 0
        %                 Pr(ii) =
        %             end
        %         end
        %%% estimate ev and svd
        %         Afun=@(u, t) ABfun(u, A, Pr, begin_in, end_in, t);
        %         [V_SVD, S_SVD, U_SVD, flag] = svds(Afun, size(A), 1, 'largest', 'Tolerance', 1e-6);
        %
        %   Compute eigenvalues of A
        %         EVAfun =@(x) B1fun(x, A, Pr, begin_in, end_in);
        %         opts.isreal = false;
        %
        %         [EV, EIG] = eigs(EVAfun, length(A), 10, 'lm', 'Tolerance', 1e-6, 'IsFunctionSymmetric', 0);
        %         [EV, ~] = qr(EV, 0);
        %%%%%%%%%%%%%%%%%%%%%%%
        rng(matInd);
        %        b = U_SVD(:,1);
        b = rand(n,1);
        %         for kk = 1 : n
        %             if(b(kk) < 0)
        %                 b(kk)=0;
        %             end
        %         end
        %         b = ones(n,1);
        b = b/norm(b);
        %         b = zeros(n, 1);
        %         b(matInd) = 1;
        %         afun = @(x)
        %         [b_U, b_sig, b_V] = svds(
        x_0 = zeros(n, 1);
        x_1 = x_0;
        r_0 = b - A * x_0;
        
        %Deflation by ND svd
%         [P_kND,Y_kND,C_kND,U_kND,resvecND] =  Deflationsvd(A,b,matrices,matInd,methodInd,methodDef,resvec,adaptive_deflation,V,H,npes,k,cycle,maxit,maxcycle,x_0,x_1,r_0,maxdef,tol,method);
        
        [P_kND,Y_kND,C_kND,U_kND,resvecND] = DeflationsvdRRQR(A,b,matInd,methodInd,methodDef,resvec,adaptive_deflation,V,H,k,cycle,maxit,maxcycle,x_0,x_1,r_0,maxdef,tol,'svdrrqr');
        
        
        %Deflation by svd
        [P_k,Y_k,C_k,U_k,resvec] = Deflationsvd(A,b,matInd,methodInd,methodDef,resvec,adaptive_deflation,V,H,k,cycle,maxit,maxcycle,x_0,x_1,r_0,maxdef,tol,method);
        
        %Deflation by EV
%         [P_kEV,Y_kEV,C_kEV,U_kEV,resvecEV] = DeflationEV(A,b,matInd,methodInd,methodDef,resvec,adaptive_deflation,V,H,k,cycle,maxit,maxcycle,x_0,x_1,r_0,maxdef,tol,'heig');
        
        
        %residual  by GMRES
        [x,flag,relres,iter,resvecGMRES] = gmres(A,b,[],tol,maxit*maxcycle,[],[],x_0);
    end
    fprintf('Method %s finished\n', methodDef{methodInd});
    plot(log10([resvecND{:, methodInd}]));
    hold on;
    plot(log10([resvec{:, methodInd}]));
%       plot(log10([resvecEV{:, methodInd}]));
    plot(log10(resvecGMRES));
    k = 0;
end
legend('Deflation based on SVD with Nested Dissection','Deflation based on SVD','GMRES')

% legend('Deflation based on SVD with Nested Dissection','Deflation based on SVD','Deflation based on eigenvectors','GMRES')
% ylim([-0.5 0])