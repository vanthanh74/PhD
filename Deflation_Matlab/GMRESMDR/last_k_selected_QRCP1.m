function last_k_columns = last_k_selected_QRCP1(A1,A2,A,k)

% ii=i
% A1 = A_prod_dissect(:,A1_global_indices_matrix{ii});
% A2 = A_prod_dissect(:,A2_global_indices_matrix{ii});
% A  =  A_prod_dissect;

A0=A;

%QRCP for first subdomain
A1=full(A1);
[Q1, R1, P1] = qr(A1,'vector');
if length(P1)>k
    [~,P1_selected_columns]=ismember(A1(:,P1(end-k+1:end)).',A0.','rows');
else
    [~,P1_selected_columns]=ismember(A1(:,P1).',A0.','rows');
end
% fprintf('Step 1 selected columns for subdomain 1')
% P1_selected_columns'
A1_last = A0(:,P1_selected_columns);

%QRCP for second subdomain
A2=full(A2);
[Q2,R2,P2] = qr(A2,'vector') ;
if length(P2)>k
    [~,P2_selected_columns]=ismember(A2(:,P2(end-k+1:end)).',A0.','rows');
else
    [~,P2_selected_columns]=ismember(A2(:,P2).',A0.','rows');
end
% fprintf('Step 1 selected columns for subdomain 2')

%remove duplicated columns 
P2_selected_columns = setdiff(P2_selected_columns,P1_selected_columns);

% P2_selected_columns'
A2_last = A0(:,P2_selected_columns);


%QRCP for B matrix
B = [A1_last A2_last];
B=full(B);
[Q,R,P_B] = qr(B,'vector');
% fprintf('Step 2 selected columns ')
[~,last_k_columns]=ismember(B(:,P_B(end-k+1:end)).',A0.','rows');
% last_k_columns'
% rejected_columns = setdiff(1:size(A,2),selected_columns_B_in_A);

%permutted A
% A_final = [A(:,rejected_columns) A(:,selected_columns_B_in_A)];
