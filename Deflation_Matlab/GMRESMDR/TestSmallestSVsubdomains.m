close all

%check if the selected columns give good approximation for the smallest
%singular values of R22
% k=20; % number of last smallest singular values
% npes=8; %number of subdomains
% A1_global_indices_matrix; % number of first subdomain
% A2_global_indices_matrix; % number of second subdomain
A_global_indices_matrix = [A1_global_indices_matrix;A2_global_indices_matrix]
% last_k_columns_matrix; % last k columns selected from npes subdomains

%%
i=8 %subdomain ith
last_k_columns_matrix_i=last_k_columns_matrix(i,:)
subdomain_i = A0(:,A_global_indices_matrix{i});
% spy(subdomain_i)

rejected_columns = setdiff(A_global_indices_matrix{i},last_k_columns_matrix_i);
subdomain_i_permuted = [A0(:,rejected_columns) A0(:,last_k_columns_matrix_i)];
% figure
% spy(subdomain_i_permuted)


% subdomain_i = A0(:,A1_global_indices_matrix{i});
% subdomain_i = A0(:,A2_global_indices_matrix{i});
[Q_i,R_i,P_i] = qr(full(subdomain_i_permuted),0);
R22_subdomain_i=R_i(end-k+1:end,end-k+1:end); % size k x k

% spy(R22_subdomain_i)

%svd of subdomain1
% ssubdomain_i=svd(full(subdomain_i));
% ssubdomain_i=svds(subdomain_i,k,'smallest');
ssubdomain_i_k=svds(subdomain_i,k,'smallest');
% ssubdomain_i_k = ssubdomain_i(end-k+1:end);

%svd of R22_subdomain1
% sR22_subdomain_i=svd(full(R22_subdomain_i));
% sR22_subdomain_i=svds(R22_subdomain_i,k,'smallest');
sR22_subdomain_i_k=svds(R22_subdomain_i,k,'smallest');
% sR22_subdomain_i_k = sR22_subdomain_i(end-k+1:end);

%compare
compare_mat = [ssubdomain_i_k sR22_subdomain_i_k]

figure
plot(1:k,ssubdomain_i_k,'db',1:k,sR22_subdomain_i_k,'.r')
legend(['subdomain',num2str(i)'],'R_{22}')
title('k smallest singular values')

%% check with original A

% %remove duplicated columns
last_k_columns_matrix_i=last_k_columns_matrix(1,:);
last_k_columns_matrix_temp = cell(number_subdomains,1);
last_k_columns_matrix_temp{1} = last_k_columns_matrix(1,:);
for i=2:2*number_subdomains
%     i=i+1
    last_k_columns_matrix_i = setdiff(last_k_columns_matrix(i,:),cat(2,last_k_columns_matrix_temp{1:i-1,:}));
    last_k_columns_matrix_temp{i} = last_k_columns_matrix_i;
end


% %select all interfaces
% all_interfaces_columns=[];
% for i=3:2+log2(npes)
% % %         i=i+1
%     for j=1:number_subdomains
% % %             j=j+1
%         if size(subdomain_indices{j,i},1)~= 0
%             all_interfaces_columns =  [all_interfaces_columns subdomain_indices(j,i)];
%         end
%     end
% end
% 
% % %remove all interface columns that already appear in selected columns in each subdomain
% all_interfaces_columns_temp=[];
% for j=1:length(all_interfaces_columns)
%     all_interfaces_columns_j = setdiff(all_interfaces_columns{j},cat(2,last_k_columns_matrix_temp{1:2*number_subdomains,:}));
%     all_interfaces_columns_temp{j} = all_interfaces_columns_j;
% end


% % select global interfaces
global_interfaces_columns =[];
for i=4:2+log2(npes)
% %     i=i+1
    for j=1:number_subdomains
% %             j=j+1
        if size(subdomain_indices{j,i},1)~= 0
%             subdomain_indices(j,i)
%             BB =  [BB A0(:,subdomain_indices{j,i})];
            global_interfaces_columns = [global_interfaces_columns subdomain_indices(j,i)];
        end
    end
end

% %remove global interface columns that already appear in selected columns in each subdomain
global_interfaces_columns_temp=[];
for j=1:length(global_interfaces_columns)
    global_interfaces_columns_j = setdiff(global_interfaces_columns{j},cat(2,last_k_columns_matrix_temp{1:2*number_subdomains,:}));
    global_interfaces_columns_temp{j} = global_interfaces_columns_j;
end


%QRCP for BB matrix
BB=[];
for i=1:2*number_subdomains
    BB = [BB A0(:,last_k_columns_matrix_temp{i})];
end



% %add all interfaces
% for i=1:length(all_interfaces_columns)
%     BB =  [BB A0(:,all_interfaces_columns_temp{i})];
% end

% %add global interfaces
for i=1:length(global_interfaces_columns)
    BB =  [BB A0(:,global_interfaces_columns_temp{i})];
end



% subdomain_indices_fullfill{end}
% figure
% spy(A0(:,subdomain_indices_fullfill{end}))
% figure
% spy(A0)
% figure
% spy(BB)

[Q,R,P_BB] = qr(full(BB),'vector');
% fprintf('Step 2 selected columns ')
[~,last_k_columns_QRCP_TournamentPivoting]=ismember(BB(:,P_BB(end-k+1:end)).',A0.','rows');

% last_k_columns_matrix_i'
% last_k_columns_matrix_ii = last_k_columns_QRCP_A;
rejected_columns = setdiff(1:size(Bcopy,2),last_k_columns_QRCP_TournamentPivoting);
A_final = [Bcopy(:,rejected_columns) Bcopy(:,last_k_columns_QRCP_TournamentPivoting)];
% P_kND = [rejected_columns last_k_columns_matrix_i'];
% figure
% spy(Bcopy)

[Q_final,R_final] = qr(full(A_final));   
% normest(A_kND-Q_final*R_final)

R22_final=R_final(end-k+1:end,end-k+1:end); % size k x k

%svd of R22_A
sR22_final   = svd(full(R22_final));
sR22_final_k = sR22_final(end-k+1:end)
% sR22_final_k=svds(R22_final,k,'smallest');
last_k_columns_QRCP_TournamentPivoting

%%



%QRCP on A
[Q_A,R_A,P_A] = qr(full(A),'vector');   

%last k columns from QRCP on A
last_k_columns_QRCP_A = P_A(end-k+1:end)
R22_A=R_A(end-k+1:end,end-k+1:end); % size k x k

%svd of R22_A
sR22_A=svd(full(R22_A));
sR22_A_k = sR22_A(end-k+1:end);


% compare_mat = [sA_k sR22_A_k]
% figure
% plot(1:k,sA_k,'db',1:k,sR22_A_k,'.r')
% legend('A','R_{22}')
% title('k smallest singular values')

compare_last_k_columns = [last_k_columns_QRCP_A; last_k_columns_QRCP_TournamentPivoting'];
compare_last_k_columns=sort(compare_last_k_columns,2)'

%---------------------------------------------------
%svd of A
sA=svd(full(A));
sA_k = sA(end-k+1:end);
% sA_k=svds(A,k,'smallest');


% sA_final=svd(full(A_final));
% sA_final_k = sA_final(end-k+1:end);

% normest(sA_k-sA_final_k)



last_k_columns_matrix_i'

%compare 
% 'SV of A',  'SV of R_{22} from QRCP of A' ,  'SV of R_{22} from QRCP with TP'
compare_mat = [ sA_k sR22_A_k sR22_final_k]

figure
plot(1:k,sA_k,'db',1:k,sR22_final_k,'.r')
legend('A','R_{22}')
title('k smallest singular values')

figure
subplot(121)
spy(R22_A)
subplot(122)
spy(R22_final)

figure
% subplot(121)
spy(A0(:,last_k_columns_QRCP_A))

% subplot(122)
figure
spy(A0(:,last_k_columns_QRCP_TournamentPivoting))
columns_QRCP_A = A0(:,last_k_columns_QRCP_A);
columns_QRCP_TP = A0(:,last_k_columns_QRCP_TournamentPivoting);
figure
spy(A0)





