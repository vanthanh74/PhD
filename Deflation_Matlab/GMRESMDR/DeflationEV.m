function [P_k,Y_k,C_k,U_k,resvec] = DeflationEV(A,b,matInd,methodInd,methodDef,resvec,adaptive_deflation,V,H,k,cycle,maxit,maxcycle,x_0,x_1,r_0,maxdef,tol,method)

if(k > 0)
            if(reduction && strcmp(matrices{matInd - 1}, matrices{matInd}) == 0)
                reduce_def_ss;
            end
%             [Q, R] = qr(A * blcLS(Pr, begin_in, end_in, Y_k), 0);
            [Q, R] = qr(A * Y_k, 0);
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
            for j = 1 : maxit
%                 w = A * blcLS(Pr, begin_in, end_in, V(:, j)); %blcLS  y(begin_in(i) : end_in(i), :) = Pr{i} \ x(begin_in(i) : end_in(i), :);
                w = A * V(:, j);
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
%                 [Q, R] = qr(A * blcLS(Pr, begin_in, end_in, Y_k), 0);
                [Q, R] = qr(A * Y_k, 0);
                C_k = Q;
                U_k = Y_k / R;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%  Singular values %%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif(strcmp(method, 'svd') == 1)
                [AP_k, Sigma_k, P_k] = svd(H,0);
                P_k = P_k(:, end - k + 1 : end);
                Y_k = V(:, 1 : maxit) * P_k;
%                 [Q, R] = qr(A * blcLS(Pr, begin_in, end_in, Y_k), 0);
                [Q, R] = qr(A * Y_k, 0);
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
%                 [Q, R] = qr(A * blcLS(Pr, begin_in, end_in, Y_k), 0);
                [Q, R] = qr(A * Y_k, 0);
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
            for j = 1 : maxit - k
%                 w = A * blcLS(Pr, begin_in, end_in, V(:, j));
                w = A * V(:, j);
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
            Vtilde = [Utilde, V(:, 1 : IndexOfLastIter)];
            Wtilde = [C_k, V(:, 1 : IndexOfLastIter + 1)];
%             G = Wtilde' * A * blcLS(Pr, begin_in, end_in, [Utilde, V(:, 1 : IndexOfLastIter)]);   
            G = Wtilde' * A * [Utilde, V(:, 1 : IndexOfLastIter)];
            [QG, RG] = qr(G);
            res = QG' * c;
            
            cycle = cycle + 1;
            y = RG(1 : end - 1, 1 : end) \ res(1 : end - 1);
            x_1 = x_1 + Vtilde * y;
%             r_1 = b - A * blcLS(Pr, begin_in, end_in, x_1);
            r_1 = b - A * x_1;
            fprintf('Real residual  = %e\n', norm(r_1));
            
            %r_1 = r_1 - Wtilde * H * y;
%             fprintf('Expanded residual  = %e\n', norm(r_1));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Harmonic Ritz values %%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            k = min(maxdef, j + k);
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
                M1 = G' * G;
                M2 = Vtilde' * Vtilde;
                [P_k, Eig] = eig(M1, M2); % M1*P_k= M2*P_k*Eig
                [~, indices] = sort(diag(Eig));
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

%         if(write_results == 1)
%             if(matInd < 10)
%                 filename = [pwd, '/results_19/', methodDef{methodInd}, '_', num2str(maxit),...
%                     '_', num2str(maxdef), '_', matrices{matInd}(11 : end-4), '_0', ...
%                     num2str(matInd), '.txt'];
%             else
%                 filename = [pwd, '/results_19/', methodDef{methodInd}, '_', num2str(maxit),...
%                     '_', num2str(maxdef), '_', matrices{matInd}(11 : end-4), '_', ...
%                     num2str(matInd), '.txt'];
%             end
%             fileid = fopen(filename, 'w+');
%             fprintf(fileid, date);
%             fprintf(fileid, '\n');
%             if(adaptive_deflation == 1)
%                 fprintf(fileid, 'Adaptive deflation: 1\n');
%             else
%                 fprintf(fileid, 'Adaptive deflation: 0\n');
%             end
%             if(reduction == 1)
%                 fprintf(fileid, 'Reduction: 1\n');
%                 if(matInd > 1 && strcmp(matrices{matInd - 1}, matrices{matInd}) == 0)
%                     fprintf(fileid, 'Tol redution = %e\n', tol_k);
%                 end
%             else
%                 fprintf(fileid, 'Reduction: 0\n');
%             end
%             fprintf(fileid, 'Iteration number: %d\n', length(resvec{matInd, methodInd}));
%             fprintf(fileid, '%E\n', resvec{matInd, methodInd}');
%             fclose(fileid);
%         end
%     end
end