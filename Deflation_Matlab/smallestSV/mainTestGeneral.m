clear
clc
close all


n = 100;
k = 10; % <n/(2*nn)

%Partition of A
nn = 2; % A= A1 U A2 U ... U Ann

kappa = 1;
[A,b,x] = heat(n,kappa);
% A= rand(n,n);

% Matlab QR with Column Pivoting
[QQ,RR,PP] = qr(A,0);

% R22_ref= Q(:,P(end-k+1:end))*R(P(end-k+1:end),:);
R22_ref = RR(end-k+1:end,end-k+1:end);
[U_ref,S_ref,V_ref] = svd(R22_ref);
[U,S,V] = svd(A);

S_diag=diag(S);
S_ref;
S_A=S_diag(end-k+1:end);
norm(S_ref-S_A,2);

S_A_last=S_A;
% QR with tournament pivoting
% A1= zeros(k,k);
% A2= zeros(k,k);


n_tp= nn/2; % number of tournament pivoting turns

A_sub={};
for i=1:n_tp
    for ii =1:nn
        % i= i+1
        A_sub{ii}= A(:,(ii-1)*n/nn+1:ii*n/nn);
    end

    %QRCP for each A_Sub
    Q={};
    R={};
    P={};
    for ii =1:nn
        [Q{ii},R{ii},P{ii}] = qr(A_sub{ii},0);
    end

    %Select the last k columns from each A_sub to make A_sub_child0
    %Sort in order each P{i}
    for ii=1:nn
        P{ii}=sort(P{ii}(n/nn-k+1:n/nn))
    end

    A_sub_child0={};
    for ii=1:nn
        A_sub_child0{ii}= A_sub{ii}(:,P{ii});
    end


    %Put the last k columns from A_sub_child0{i} together
     A_sub_child1={};
    for ii=1:n_tp
        A_sub_child1=[A_sub_child0{(ii-1)*n_tp+1} A_sub_child0{(ii-1)*n_tp+2}]
    end
    A_sub=A_sub_child1;
    nn=nn/2;
end


% Q2={};
% R2={};
% P2={};
% for i=1:n_tp
%      [Q2{i},R2{i},P2{i}] = qr(A_sub_child1{i},0);
% end
% 
% for i=1:n_tp
%     P2{i}=sort(P2{i}(n/nn-k+1:n/nn))
% end

% A_sub_child2={};

% A_sub_child2= A_sub(:,P{i})


% [Q0,R0,P0] = qr(A_sub_child1,0);
% P0=sort(P0(end-k+1:end));
% A0 =  A_sub_child1(:,P0);

[U_TP,S_TP,V_TP] = svd(A_sub_child1);
S_TP_dia=diag(S_TP);
S_TP_last = S_TP_dia(end-k+1:end);
S_A;
S_A_last
S_TP_last
norm(S_A_last-S_TP_last,2)

figure
semilogy(1:k,S_A_last,'r',1:k,S_TP_last,'b')
legend('Singular Values of A','Singular Values of A_3')