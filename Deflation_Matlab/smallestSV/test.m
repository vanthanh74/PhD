clear
clc
close all


n = 200;
k = 5; % <n/(2*nn)

%Partition of A
nn = 3; % A= A1 U A2 U ... U Ann

kappa = 1;
% [A,b,x] = gravity(n);
[A,b,x] = heat(n,kappa);
% [A,b,x] = shaw(n); % works
% [A,b,x] = spikes(n); % works
% [A,b] = ursell(n);
% [A,b,x] = wing(n);
% [A,b,x,t] = i_laplace(n);% works
% [A,b,x] = phillips(n);
% [A,b] = parallax(n);
% [A,b,x] = baart(n);
% [A,b,x] = deriv2(n);
% [A,b,x] = foxgood(n);  % works




% A= rand(n, n);
% rank(A,10^-3)
% Matlab QR with Column Pivoting
[QQ,RR] = qr(A);

% R22_ref= Q(:,P(end-k+1:end))*R(P(end-k+1:end),:);
R22_ref = RR(end-k+1:end,end-k+1:end); 
% [U_ref,S_ref,V_ref] = svd(R22_ref);
[U,S,V] = svd(A);
S_A_last = diag(S(end-k+1:end,end-k+1:end));
diag_R22_ref = diag(R22_ref);         % approximate for the last 
                                       %k smallest singular values of A
% S_diag=diag(S);
% S_ref;
% S_A=S_diag(end-k+1:end);
norm(S_A_last-diag_R22_ref,2)

% S_A_last=S_A;
% S_ref_diag = diag(S_ref);
% S_ref_last = S_ref_diag(end-k+1:end);
% QR with tournament pivoting
% A1= zeros(k,k);
% A2= zeros(k,k);


n_tp= nn-1; % number of tournament pivoting turns
nn0 = nn;
A_sub={};
for ii =1:nn  
    A_sub{ii}= A(:,(ii-1)*n/nn+1:ii*n/nn);
end

for i=1:n_tp  
    % ii= ii+1
    % i= i+1
    %QRCP for each A_Sub
    Q={};
    R={};
    P={};
    for ii =1:nn
        [Q{ii},R{ii},P{ii}] = qr(A_sub{ii},0);
    end

    
    % determine the last k columns from each P{i}
    for ii =1:nn
        P{ii} = P{ii}((end-k+1:end));
    end
    
    
    %Sort in order each P{i}
%      P{ii}=sort( P{ii})
%     for ii=1:nn
%         P{ii}=sort(P{ii}(length(P{ii})-k+1:length(P{ii})));
%     end
    
    %Select the last k columns from each A_sub to make A_sub_child0
    A_sub_child0={};
    for ii=1:nn
        A_sub_child0{ii}= A_sub{ii}(:,P{ii});
    end


    %Put the last k columns from A_sub_child0{i} together
     A_sub_child1={};
    for ii=1:n_tp
        % ii = ii+1
        A_sub_child1{ii}=[A_sub_child0{(ii-1)*2+1} A_sub_child0{(ii-1)*2+2}];
    end
    A_sub=A_sub_child1;
    nn=nn/2;
   
    if n_tp==1
        break
    end
     n_tp=nn/2;
end



A_sub_child1=A_sub_child1{1};

[Q0,R0,P0] = qr(A_sub_child1,0);
P00=sort(P0(end-k+1:end));
% A_sub_child1(:,P00)
% [R,C]=find(A==A_sub_child1(:,1))
[~,P00_in_A]=ismember(A_sub_child1(:,P00).',A.','rows'); % find the corresponding worst columns in A
P00_in_A=P00_in_A'
% P00_in_A=sort(P00_in_A)'
% find(ismember(A_sub_child1(:,1),A),n)
% P01=sort(P0(1:end-k));
% A0 =  A_sub_child1(:,P0);

P01_in_A = 1:n;
P01_in_A(:,P00_in_A)=[];

% [Q3,R3,P3] = qr(A3,0);
% P31=sort(P3(end-k+1:end));
% P30=sort(P3(1:end-k));
% A33 =  A3(:,P31);

% Put the worst columns of A to the end
A_permuted= [A(:,P01_in_A) A(:,P00_in_A)];  % [A(:,P01) A(:,P00)]

[Q4,R4] = qr(A_permuted);
% [Q44,R44,P44] = qr(A,0);

R22_TP = R4(end-k+1:end,end-k+1:end);
% R22_44TP = R44(end-k+1:end,end-k+1:end);

% R22_PT = A3(end-k+1:end,end-k+1:end)

% [U_TP,S_TP,V_TP] = svd(R22_TP);



% [U_TP,S_TP,V_TP] = svd(A_sub_child1);
S_TP_diag = diag(R22_TP);
S_TP_last = S_TP_diag(end-k+1:end);
% S_A;
% S_A_last
S_TP_last
S_A_last
diag_R22_ref
% S_ref_last
% norm(S_A_last-S_TP_last,2)


figure
plot(1:k,diag_R22_ref,'-sr',1:k,S_TP_last,'-^b',1:k,S_A_last,'-dm')
legend('Last smallest Singular Values of R_{22}ref-QR(A)','Last smallest Singular Values of A_{permuted}','Last smallest Singular Values of A-SVD(A)')

figure
plot(1:k,diag_R22_ref,'-sr',1:k,S_TP_last,'-^b')
legend('Last smallest Singular Values of R_{22}ref-QR(A)','Last smallest Singular Values of A_{permuted}')


figure
plot(1:k,S_TP_last,'-^b',1:k,S_A_last,'-dm')
legend('Last smallest Singular Values of A_{permuted}','Last smallest Singular Values of A-SVD(A)')
