clear all; clc; close all;
format short e;

% save results
write_results = 0;

% addpath('matrix');
% The test matrices
% Pr_BIGP0='matrix/BIGP1_METRIC1_ELX_TSTEP_114_NEWTONIT_0_Pressure.mtx';
% Pr_BIGP1='matrix/BIGP1_METRIC1_ELX_TSTEP_114_NEWTONIT_1_Pressure.mtx';

% Pr_BIGP0='fidap001.mtx';%N=317
Pr_BIGP0='fidap003.mtx';%N=1821
% Pr_BIGP0='fidap021.mtx';%N=656
% Pr_BIGP0='fidap024.mtx';%N=2283 nz=47897
% Pr_BIGP0='g7jac040sc.mtx';
% Pr_BIGP0='cavity01.mtx';
% Pr_BIGP0='dw4096.mtx';%N=8192 nz=41746

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
[p1, edgecut, s, begin_in, end_in] = DefineAlpha(A1, nb_subdomains);
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
A = A1(p1,p1);
n = length(A1);


figure
subplot(121)
spy(A1(p1,p1))
subplot(122)
spy(A1(p1,p1))

figure
spy(A)





% figure
% spy(A(begin_in(1):end_in(1),begin_in(1):end_in(1)))




% Block Jacobi
Pr = cell(nb_subdomains, 1);

% Parameters of GMRES
maxit = 25;%50;%
maxdef = 5;%10;%
maxcycle = 100;
k = 0;
rng(1);
% Matrices of Krylov
H = zeros(maxit + 1, maxit);
V = zeros(n, maxit + 1);
resvec = cell(length(matrices), 3);

% Method of recycling
adaptive_deflation = 0;
reduction = 0;
methodDef = {'svd','reig', 'heig'};
      methodDef = {'svd'};
for methodInd = 1 : length(methodDef)
    method = methodDef{methodInd};
    
    for matInd = 1 : length(matrices)
        cycle = 1;
        if (matInd == 1)
            tol = tol_1;
            if( adaptive_deflation == 1 )
               method = 'reig';
            end
        else
            tol = tol_2;
            method = methodDef{methodInd};
        end
        resvec{matInd, methodInd} = [];
        A1 = mmread(matrices{matInd});
        [A1,R,C] = equil_rar(A1);
        
        A = A1(p1, p1);
%         figure
%         spy(A)
        for i = 1 : length(begin_in)
            if begin_in(i)*end_in(i)~=0
                Pr{i} = A(begin_in(i) : end_in(i), begin_in(i) : end_in(i));
            end
        end
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
        if(k > 0)
            if(reduction && strcmp(matrices{matInd - 1}, matrices{matInd}) == 0)
                reduce_def_ss;
            end
            [Q, R] = qr(A * blcLS(Pr, begin_in, end_in, Y_k), 0);
            C_k = Q;
            U_k = Y_k /R;
            x_1 = x_0 + U_k * (C_k' * r_0);
            r_1 = r_0 - C_k * (C_k' * r_0);
            norm(C_k' * r_0)
        else
            normr = norm(r_0);
            V(:, 1) = r_0/normr;
            e_1 = zeros(maxit + 1, 1);
            e_1(1) = 1;
            c = normr * e_1;
            for j = 1 : maxit  %perform maxit steps of gmres, solving min||c-H_my||_2 for y 
                w = A * blcLS(Pr, begin_in, end_in, V(:, j)); % and generating V_{m+1} and H_m
                temp = V(:, 1 : j)' * w;
                w = w - V(:, 1 : j) * temp;
                temp1 = V(:, 1 : j)' * w;
                w = w - V(:, 1 : j) * temp1;
                H(1 : j, j) = temp + temp1;
                H(j + 1, j) = norm(w);
                V(:, j + 1) = w / H(j + 1, j);
                [QH, RH] = qr(H(1 : j + 1, 1 : j));
                e_1 = zeros(j + 1, 1);
                e_1(1) = 1;
                c = normr * e_1;
                res = QH' * c;
                fprintf('Res(%d) = %e\n', j, abs(res(end)));
                resvec{matInd, methodInd} = [resvec{matInd, methodInd}, abs(res(end))];
            end
            if( adaptive_deflation == 1 )
                if( abs(res(end))/normr > 8/9)
                    method = 'reig';
                    fprintf('method adaptively changed from %s to reig\n', methodDef{methodInd});
                else
                    method = methodDef{methodInd};
                end
            end
            k = maxdef;
            y = RH(1 : end - 1, :) \ res(1 : end - 1);
            x_1 = x_0 + V(:, 1 : maxit) * y;
            r_1 = V * (c - H * y);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Harmonic Ritz values %%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if(strcmp(method, 'heig') == 1)
                Hm = H(1 : maxit, 1 : maxit);
                hm = H(maxit + 1, maxit);
                em = zeros(maxit, 1);
                em(maxit) = 1;
                E = Hm + hm * hm * (Hm \ em * em');
                [P_k, Eig] = eig(E);
                [~, indices] = sort(diag(Eig));
                P_k = P_k(:, indices(1 : k));
                Y_k = V(:, 1 : maxit) * P_k;
                [Q, R] = qr(A * blcLS(Pr, begin_in, end_in, Y_k), 0);
                C_k = Q;
                U_k = Y_k / R;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%  Singular values %%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif(strcmp(method, 'svd') == 1)
                [AP_k, Sigma_k, P_k] = svd(H,0); % compute k deflation vectors P_k
%                 size(H) size(A1)
                 % compute k deflation vectors P_k using Nested Dissection QRCP 
                 % [U_k,Sigma_k,V_k] = NDQRCPSingularValuesApproximation(A1,npes,sizes,p1,k)
                [AP_k,Sigma_k,P_k,V] = NDQRCPSingularValuesApproximation(A1,b,Pr,begin_in,end_in,k,maxit);
%                 P_k1 = P_k1(:, end - k + 1 : end);
%                 return
                P_k = P_k(:, end - k + 1 : end);

                Y_k = V(:, 1 : maxit) * P_k;
                [Q, R] = qr(A * blcLS(Pr, begin_in, end_in, Y_k), 0);
                C_k = Q;
                U_k = Y_k / R;
                %         Y_k = V(:, 1 : maxit) * P_k(:, end - k + 1: end);
                %         C_k = V * AP_k(:, end - k + 1: end);
                %         U_k = Y_k / Sigma_k(end - k : end - 1, end - k + 1 : end);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%  Ritz values %%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif(strcmp(method, 'reig') == 1)
                Hm = H(1 : maxit, 1 : maxit);
                [P_k, Eig] = eig(Hm);
                [~, indices] = sort(diag(Eig));
                P_k = P_k(:, indices(1 : k));
                Y_k = V(:, 1 : maxit) * P_k;
                [Q, R] = qr(A * blcLS(Pr, begin_in, end_in, Y_k), 0);
                C_k = Q;
                U_k = Y_k / R;
            end
        end
        residual = 1;
        while(cycle < maxcycle && norm(r_1) > tol)
%             fprintf('Cycle number %d\n', cycle);
%             method = methodDef{mod(cycle,3) + 1};
            normr = norm(r_1);
            H = 0 * H;
            V = 0 * V;
            ctr = C_k' * r_1;
            r_1 = r_1 - C_k * ctr;
            x_1 = x_1 + U_k * ctr;
            V(:, 1) = r_1/normr;
            D = zeros(k, k);
            %D = diag(diag(Sigma_k(end - k : end - 1, end - k + 1 : end)));
            for i = 1 : k
                D(i, i) = 1/norm(U_k(:, i));
            end
            Utilde = U_k * D;
            H(1 : k, 1 : k) = D;
%             return
            for j = 1 : maxit - k
%                 j=j+1
                w = A * blcLS(Pr, begin_in, end_in, V(:, j));
                temp = C_k' * w;
                w = w - C_k * temp;
                temp2 = C_k' * w;               % double ...
                w = w - C_k * temp2;            % orthogonalization
                H(1 : k, k + j) = temp + temp2;
                temp = V(:, 1 : j)' * w;
                w = w - V(:, 1 : j) * temp;
                temp1 = V(:, 1 : j)' * w;
                w = w - V(:, 1 : j) * temp1;
                H(k + 1 : k + j, k + j) = temp + temp1;
                H(k + j + 1, k + j) = norm(w);
                V(:, j + 1) = w / H(k + j + 1, k + j);
                [QH, RH] = qr(H(1 : k + j + 1, 1 : k + j));
                c = [C_k V(:, 1 : j + 1)]' * r_1;
                res = QH' * c;
                fprintf('Res(%d) = %e\n', j, abs(res(end)));
                residual = abs(res(end));
                resvec{matInd, methodInd} = [resvec{matInd, methodInd}, abs(res(end))];
                IndexOfLastIter = j;
                if(residual < tol)
                    break;
                end
            end
            if( adaptive_deflation )
                if( residual/normr > 2/3 && j >= maxit - k)
                    method = 'reig';
                    if(residual/normr > 0.99 && strcmp( methodDef{methodInd}, 'reig') == 1 )
                        method = 'heig';
                    end
                    fprintf('method adaptively changed from %s to %s\n', methodDef{methodInd}, method);
                else
                    method = methodDef{methodInd};
                end
            end
%             return
            Vtilde = [Utilde, V(:, 1 : IndexOfLastIter)];
            Wtilde = [C_k, V(:, 1 : IndexOfLastIter + 1)];
            G = Wtilde' * A * blcLS(Pr, begin_in, end_in, [Utilde, V(:, 1 : IndexOfLastIter)]);
            [QG, RG] = qr(G);
            res = QG' * c;
            
            cycle = cycle + 1;
            y = RG(1 : end - 1, 1 : end) \ res(1 : end - 1);
            x_1 = x_1 + Vtilde * y;
            r_1 = b - A * blcLS(Pr, begin_in, end_in, x_1);
            fprintf('Real residual  = %e\n', norm(r_1));
            
            %r_1 = r_1 - Wtilde * H * y;
%             fprintf('Expanded residual  = %e\n', norm(r_1));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Harmonic Ritz values %%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            k = min(maxdef, j + k);
%             return
            if(strcmp(method, 'heig') == 1)
                M1 = G' * G;
                M2 = G' * (Wtilde' * Vtilde);
                [P_k, Eig] = eig(M1, M2);
                [~, indices] = sort(diag(Eig));
                P_k = P_k(:, indices(1 : k));
                Y_k = Vtilde * P_k;
                [Q, R] = qr(G * P_k, 0);
                C_k = Wtilde * Q;
                U_k = Y_k / R;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%  Singular values %%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif(strcmp(method, 'svd') == 1)
                M1 = G' * G;              %   [V,D] = EIG(A,B) produces a diagonal matrix D of generalized                                                                     
                M2 = Vtilde' * Vtilde;    %   eigenvalues and a full matrix V whose columns are the corresponding
%                 [P_k, Eig] = eig(M1, M2); %   eigenvectors so that A*V = B*V*D.                                   
%                 [~, indices] = sort(diag(Eig));
                
                [AP_k,Sigma_k,P_k,V] = NDQRCPSingularValuesApproximationwithJacobiPr(A1,b,Pr,begin_in,end_in,npes,k,maxit);                 
                [~, indices] = sort(diag(Sigma_k));
%                  P_k = P_k(:, end - k + 1 : end);
%                            size(Vtilde)
%                            size(P_k)
                P_k = P_k(:, indices(1 : k));
                Y_k = Vtilde * P_k;
                [Q, R] = qr(G * P_k, 0);
                C_k = Wtilde * Q;
                U_k = Y_k / R;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%  Ritz values %%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif(strcmp(method, 'reig') == 1)
                M1 = Vtilde' * Wtilde * G;
                M2 = Vtilde' * Vtilde;
                [P_k, Eig] = eig(M1, M2);
                [~, indices] = sort(diag(Eig));
                P_k = P_k(:, indices(1 : k));
                Y_k = Vtilde * P_k;
                [Q, R] = qr(G * P_k, 0);
                C_k = Wtilde * Q;
                U_k = Y_k / R;
            end
            %%%%%%%% Harmonic ritz %%%%%%%%%%%
        end
        fprintf('Method %s system %d converged in %d iteration done over %d cycles\n',...
            methodDef{methodInd}, matInd, length(resvec{matInd, methodInd}), cycle);

%         plot(log10(resvec{matInd}));
%         hold on;
        if(write_results == 1)
            if(matInd < 10)
                filename = [pwd, '/results_19/', methodDef{methodInd}, '_', num2str(maxit),...
                    '_', num2str(maxdef), '_', matrices{matInd}(11 : end-4), '_0', ...
                    num2str(matInd), '.txt'];
            else
                filename = [pwd, '/results_19/', methodDef{methodInd}, '_', num2str(maxit),...
                    '_', num2str(maxdef), '_', matrices{matInd}(11 : end-4), '_', ...
                    num2str(matInd), '.txt'];
            end
            fileid = fopen(filename, 'w+');
            fprintf(fileid, date);
            fprintf(fileid, '\n');
            if(adaptive_deflation == 1)
                fprintf(fileid, 'Adaptive deflation: 1\n');
            else
                fprintf(fileid, 'Adaptive deflation: 0\n');
            end
            if(reduction == 1)
                fprintf(fileid, 'Reduction: 1\n');
                if(matInd > 1 && strcmp(matrices{matInd - 1}, matrices{matInd}) == 0)
                    fprintf(fileid, 'Tol redution = %e\n', tol_k);
                end
            else
                fprintf(fileid, 'Reduction: 0\n');
            end
            fprintf(fileid, 'Iteration number: %d\n', length(resvec{matInd, methodInd}));
            fprintf(fileid, '%E\n', resvec{matInd, methodInd}');
            fclose(fileid);
        end
    end
    fprintf('Method %s finished\n', methodDef{methodInd});
    plot(log10([resvec{:, methodInd}]));
    hold on;
    k = 0;
end