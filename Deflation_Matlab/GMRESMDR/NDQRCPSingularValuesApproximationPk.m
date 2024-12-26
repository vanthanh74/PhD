
function [A_ND,A_permutted,P_kND,p,last_k_columns_matrix,last_k_columns_matrix_temp] = NDQRCPSingularValuesApproximationPk(A,npes,k,k_df)

% A = delsq(numgrid('L',502)); 
% k = 10; %number of selected columns


% filename = 'fidap001.mtx'%;d1=63;d2=161;d3=216;% N=216,Nz=4339,A1=1:63,A2=64:161,A3=162:216
% p=[160,159,143,142,191,161,189,188,190,195,196,194,193,192,197,144,145,140,138,137,141,139,173,153,186,180,185,179,182,187,181,214,210,212,213,215,211,216,208,209,206,207,183,152,151,184,178,177,176,154,172,170,205,201,202,203,204,171,175,200,199,198,174,69,70,5,68,98,97,4,23,25,1,3,2,66,65,64,67,35,34,24,33,32,31,30,22,29,19,93,94,92,62,63,57,58,59,16,9,85,86,10,8,26,6,11,7,21,20,28,17,18,13,14,27,12,15,61,41,43,53,52,55,54,60,56,49,48,42,39,40,44,37,36,38,51,45,46,47,50,87,95,96,88,83,89,84,91,90,79,76,78,77,72,71,74,75,73,80,81,82,131,101,99,100,165,168,169,167,166,130,118,132,109,163,162,164,129,116,117,115,103,119,158,157,156,155,110,112,108,107,106,105,104,147,146,150,149,148,120,121,102,126,113,127,114,128,123,124,125,122,135,134,136,133,111];

% filename = 'cavity01.mtx';%N=317

% filename = 'e05r0000.mtx'; %N=236

% filename = 'fidapm02.mtx'; %N=537

% filename = 'fidap001.mtx'; %N=216

% filename = 'fidap022.mtx'; %N=839

% filename = 'fidap021.mtx'; %N=656

% filename = 'fidap003.mtx'; %N=1821

% filename = 'fidap024.mtx'; %N=2283 nz=47897

% filename = 'lshp3466.mtx'; %N=3466 nz=23896

% filename = 'dw4096.mtx'; %N=8192 nz=41746

% filename = 'poli_large.mtx'; %N=13935 nz=63307

% filename = 'circuit_3.mtx'; %N=12127 nz=48137

% filename = 'g7jac040sc.mtx'; %N=11790 nz=107383

% filename = '3dtube.mtx'; %N=45330 nz=3213618

% filename = 'venkat50.mtx'; %N=62424 nz=

% filename = 'k3plates.mtx'; %N=11107 nz=378927

% filename = 'pkustk07.mtx'; %N=16,860 nz=2,418,804 

% load('bcspwr10.mat');
% A = Problem.A;
% spy(A)
% A = mmread(filename);
% fprintf('Size of A - %d\n', length(A));
% fprintf('Non zeros  of A - %d\n', nnz(A));

% type = 'BinaryTree';%
type = 'kN';%
n= length(A);
% A_prod = A'*A;
% number_iter = 5;
% npes=8;
% i=1:sqrt(npes/2);
% npes_set=2.^i;

% k_smallest_sv_approx_matrix = zeros(number_iter,k);
% apply nested dissection for A'A to get the permutation vector p
% * [perm,iperm,sizes] = metismex('NodeNDP',A,npes,options)
% for iter=1:number_iter
% iter=1
% npes_vect(iter)=npes;
[p,ip,sizes] = metismex('NodeNDP',A,npes); % 
% sizes % sizes of separators
% [t,q] = etree(A_prod(p,p));
%  p = p(q);
% figure
% spy(A(p,p))
% figure
% spy(A(ip,ip))
% npes = 32
%      
% sizes =
%   Columns 1 through 13
%     88    87    72    86    83    83    92    78    80    89    83    83    97
%   Columns 14 through 26
%     83    71    78    82    82    82    81    86    95    78    78    83    84
%   Columns 27 through 39
%     79    80    83    83    74    91    21    19    20    16    20    23    17
%   Columns 40 through 52
%     12    18    19    15    18    15    23    14    16    25    26    24    35
%   Columns 53 through 63
%     32    26    29    33    32    33    34    38    87    93    77
%    
%    
%    subdomains_sizes_matrix =
%     88    87    21     0     0     0     0
%     72    86    19    25     0     0     0
%     83    83    20     0     0     0     0
%     92    78    16    26    32     0     0
%     80    89    20     0     0     0     0
%     83    83    23    24     0     0     0
%     97    83    17     0     0     0     0
%     71    78    12    35    33    87     0
%     82    82    18     0     0     0     0
%     82    81    19    32     0     0     0
%     86    95    15     0     0     0     0
%     78    78    18    26    34     0     0
%     83    84    15     0     0     0     0
%     79    80    23    29     0     0     0
%     83    83    14     0     0     0     0
%     74    91    16    33    38    93    77
   
   
A_prod_dissect = A(p,p);

% I=eye(n,n);
% P_mat = I(:,p);
% A_prod_dissect1 = P_mat'*A*P_mat;

%  figure
% spy(A_prod_dissect1)
% title('A with ND of A^TA')
% t=tree
%%
number_subdomains=npes/2;%((2*npes-1)-1)/3;
% [p1,ip1,sizes1] = metismex('NodeNDP',A_prod,3);
subdomains_sizes_matrix = zeros(number_subdomains,2+log2(npes)); % 2+n   %(npes=2^n)
ii=0;
for i=1:number_subdomains 
    subdomains_sizes_matrix(i,1:2) = [sizes(i+ii) sizes(i+1+ii)];% sizes(2*number_subdomains+i)];
    ii=ii+1;
%     i=i+1
end
% sizes_matrix

% separators = cell(log2(npes),1)
% separators{1} = sizes(npes+1:npes+number_subdomains)
temp_no=number_subdomains;
kk=0;
count=0;
for i= 1:log2(npes)%npes+1:2*npes-1
    % i=i+1
%     if i<3*npes/2+1
%         sizes_matrix(1,log2(npes)+kk)
        separators{i} = sizes(npes+1+kk:npes+temp_no+kk);
%         separators{i}
        if length(separators{i})< number_subdomains && length(separators{i})>1            
%             subdomains_sizes_matrix(2*(count+1):2*(count+1):size(subdomains_sizes_matrix,1),rank(subdomains_sizes_matrix)+i-count-1)= separators{i}';           
            subdomains_sizes_matrix(npes/size(separators{i-1},2):npes/size(separators{i-1},2):size(subdomains_sizes_matrix,1),sum( all( subdomains_sizes_matrix ~= 0, 1)) +1)= separators{i}';                      
%             count=count+1;
        elseif length(separators{i})==1
            subdomains_sizes_matrix(end,end)= separators{i}';
        else
            subdomains_sizes_matrix(:,2+i)= separators{i}';
        end
        kk=kk+temp_no;
        temp_no=temp_no/2;
%     end
    
end

% subdomains_sizes_matrix






 %%
%  [t1,q1] = etree(A_prod(p1,p1));
%  p1 = p1(q1);
 
%  figure
% % subplot(121)
% % spy(A)
% % title('Original A')
% % % figure
% % subplot(122)
% spy(A_prod_dissect(p,p))
% title('A with ND of A^TA')

% handle subdomains
subdomain_indices = cell(number_subdomains,2+log2(npes));%zeros(number_subdomains,3);
xx=1;
for i=1:number_subdomains
    for j=1:2+log2(npes)
        if subdomains_sizes_matrix(i,j)~=0
                subdomain_indices{i,j} = xx:xx+subdomains_sizes_matrix(i,j)-1;
                xx= xx + subdomains_sizes_matrix(i,j);
        end
    end
    
end
% subdomain_indices
%%
% n=3;
% figure
% spy(A_prod_dissect(subdomain_indices{n,1}:subdomain_indices{4,5}(end),subdomain_indices{n,1}:subdomain_indices{4,5}(end)))
% % spy(A_prod_dissect(subdomain_indices{2,:},subdomain_indices{2,:}))
%%
% % subdomains grouping
% subdomains_grouping_matrix = zeros(npes,2+log2(npes)-1);
% 
% first_subdomain_local_matrices = subdomains_sizes_matrix(:,[1 3:2+log2(npes)]);
% second_subdomain_local_matrices = subdomains_sizes_matrix(:,[2 3:2+log2(npes)]);
% %fill-in interfaces
% for i=1:number_subdomains-1
% %     i=i+1
%     for j=1:1+log2(npes)
% %         j=j+1
%         if first_subdomain_local_matrices(i,j)==0 && first_subdomain_local_matrices(i+1,j)~=0
%             first_subdomain_local_matrices(i,j)= first_subdomain_local_matrices(i+1,j);
%         end
%          if second_subdomain_local_matrices(i,j)==0 && second_subdomain_local_matrices(i+1,j)~=0
%             second_subdomain_local_matrices(i,j)= second_subdomain_local_matrices(i+1,j);
%          end
%          if j==1+log2(npes)
%             first_subdomain_local_matrices(i,j) = first_subdomain_local_matrices(end,j);
%             second_subdomain_local_matrices(i,j) = second_subdomain_local_matrices(end,j);
%         end      
%     end
% end
% first_subdomain_local_matrices
% second_subdomain_local_matrices

%fill-in subdomain indices with interfaces
subdomain_indices_fullfill  = subdomain_indices;%cell(number_subdomains,2+log2(npes))

if length(separators)>1
    temp_no=number_subdomains/length(separators{2});
else
    temp_no=number_subdomains;
end
kk=1;
% ii=2;
separator_jump=2;
% find(first_subdomain_local_matrices(1,:)==subdomains_sizes_matrix(1,:))
for j=4:2+log2(npes)
%     j=j+1
    for i=1:number_subdomains
%     i=i+1          
        if size(subdomain_indices_fullfill{i,j},1)==1 && size(subdomain_indices_fullfill{i-1,j},1)==0 && mod(i,2)==0
             subdomain_indices_fullfill(kk:i-1,j)= subdomain_indices_fullfill(i,j);       
%              ii=npes_set(i);
        end
        if mod(i,separator_jump)==0
            kk=kk+temp_no;
        end
        if j==2+log2(npes)
            subdomain_indices_fullfill(i,j) = subdomain_indices_fullfill(end,j);
        end         
    end
    kk=1;
    temp_no=temp_no*2;
    separator_jump=separator_jump*2;
end
% subdomain_indices_fullfill

first_subdomain_local_matrices_indices = subdomain_indices_fullfill(:,[1 3:end]);
second_subdomain_local_matrices_indices = subdomain_indices_fullfill(:,[2 3:end]);

% first_subdomain_local_matrices_indices =
%   16×6 cell array
%     [1×88 double]    [1×21 double]    [1×25 double]    [1×32 double]    [1×87 double]    [1×77 double]
%     [1×72 double]    [1×19 double]    [1×25 double]    [1×32 double]    [1×87 double]    [1×77 double]
%     [1×83 double]    [1×20 double]    [1×26 double]    [1×32 double]    [1×87 double]    [1×77 double]
%     [1×92 double]    [1×16 double]    [1×26 double]    [1×32 double]    [1×87 double]    [1×77 double]
%     [1×80 double]    [1×20 double]    [1×24 double]    [1×33 double]    [1×87 double]    [1×77 double]
%     [1×83 double]    [1×23 double]    [1×24 double]    [1×33 double]    [1×87 double]    [1×77 double]
%     [1×97 double]    [1×17 double]    [1×35 double]    [1×33 double]    [1×87 double]    [1×77 double]
%     [1×71 double]    [1×12 double]    [1×35 double]    [1×33 double]    [1×87 double]    [1×77 double]
%     [1×82 double]    [1×18 double]    [1×32 double]    [1×34 double]    [1×93 double]    [1×77 double]
%     [1×82 double]    [1×19 double]    [1×32 double]    [1×34 double]    [1×93 double]    [1×77 double]
%     [1×86 double]    [1×15 double]    [1×26 double]    [1×34 double]    [1×93 double]    [1×77 double]
%     [1×78 double]    [1×18 double]    [1×26 double]    [1×34 double]    [1×93 double]    [1×77 double]
%     [1×83 double]    [1×15 double]    [1×29 double]    [1×38 double]    [1×93 double]    [1×77 double]
%     [1×79 double]    [1×23 double]    [1×29 double]    [1×38 double]    [1×93 double]    [1×77 double]
%     [1×83 double]    [1×14 double]    [1×33 double]    [1×38 double]    [1×93 double]    [1×77 double]
%     [1×74 double]    [1×16 double]    [1×33 double]    [1×38 double]    [1×93 double]    [1×77 double]
%     
% second_subdomain_local_matrices_indices =
%   16×6 cell array
%     [1×87 double]    [1×21 double]    [1×25 double]    [1×32 double]    [1×87 double]    [1×77 double]
%     [1×86 double]    [1×19 double]    [1×25 double]    [1×32 double]    [1×87 double]    [1×77 double]
%     [1×83 double]    [1×20 double]    [1×26 double]    [1×32 double]    [1×87 double]    [1×77 double]
%     [1×78 double]    [1×16 double]    [1×26 double]    [1×32 double]    [1×87 double]    [1×77 double]
%     [1×89 double]    [1×20 double]    [1×24 double]    [1×33 double]    [1×87 double]    [1×77 double]
%     [1×83 double]    [1×23 double]    [1×24 double]    [1×33 double]    [1×87 double]    [1×77 double]
%     [1×83 double]    [1×17 double]    [1×35 double]    [1×33 double]    [1×87 double]    [1×77 double]
%     [1×78 double]    [1×12 double]    [1×35 double]    [1×33 double]    [1×87 double]    [1×77 double]
%     [1×82 double]    [1×18 double]    [1×32 double]    [1×34 double]    [1×93 double]    [1×77 double]
%     [1×81 double]    [1×19 double]    [1×32 double]    [1×34 double]    [1×93 double]    [1×77 double]
%     [1×95 double]    [1×15 double]    [1×26 double]    [1×34 double]    [1×93 double]    [1×77 double]
%     [1×78 double]    [1×18 double]    [1×26 double]    [1×34 double]    [1×93 double]    [1×77 double]
%     [1×84 double]    [1×15 double]    [1×29 double]    [1×38 double]    [1×93 double]    [1×77 double]
%     [1×80 double]    [1×23 double]    [1×29 double]    [1×38 double]    [1×93 double]    [1×77 double]
%     [1×83 double]    [1×14 double]    [1×33 double]    [1×38 double]    [1×93 double]    [1×77 double]
%     [1×91 double]    [1×16 double]    [1×33 double]    [1×38 double]    [1×93 double]    [1×77 double]

%%
% n=1;
% % A_prod_dissect=full(A_prod_dissect);
% % figure
% % spy(A_prod_dissect([first_subdomain_local_matrices_indices{1,1} first_subdomain_local_matrices_indices{1,2} first_subdomain_local_matrices_indices{1,3} first_subdomain_local_matrices_indices{1,4}],[first_subdomain_local_matrices_indices{1,1} first_subdomain_local_matrices_indices{1,2} first_subdomain_local_matrices_indices{1,3} first_subdomain_local_matrices_indices{1,4}]))
% % figure
% % spy(A_prod_dissect([second_subdomain_local_matrices_indices{1,1} second_subdomain_local_matrices_indices{1,2} second_subdomain_local_matrices_indices{1,3} second_subdomain_local_matrices_indices{1,4}],[second_subdomain_local_matrices_indices{1,1} second_subdomain_local_matrices_indices{1,2} second_subdomain_local_matrices_indices{1,3} second_subdomain_local_matrices_indices{1,4}]))
% figure
% spy(A_prod_dissect(:,cat(2,first_subdomain_local_matrices_indices{n,:})))
% figure
% spy(A_prod_dissect(:,cat(2,second_subdomain_local_matrices_indices{n,:})))

%%
%Nested dissection for A11'A11
A1_global_indices_matrix=cell(number_subdomains,1);
for i=1:number_subdomains
    A1_global_indices_matrix{i,:} = cat(2,first_subdomain_local_matrices_indices{i,:});
end



%Nested dissection for A22'A22
A2_global_indices_matrix=cell(number_subdomains,1);
for i=1:number_subdomains
    A2_global_indices_matrix{i,:} = cat(2,second_subdomain_local_matrices_indices{i,:});
end

%combination of interfaces of subdomains
% %combination of interfaces of first subdomain
% interface_global_first_subdomain = cat(2,subdomain_indices_fullfill{1,3:end});
% figure
% spy(A_prod_dissect(:,interface_global_first_subdomain))
% 
% %combination of interfaces of second subdomain
% interface_global_second_subdomain = cat(2,subdomain_indices_fullfill{2,3:end});
% figure
% spy(A_prod_dissect(:,interface_global_first_subdomain))
interface_global_subdomain=cell(number_subdomains,1);
for i=1:number_subdomains
    interface_global_subdomain{i} = cat(2,subdomain_indices_fullfill{i,3:end});
end

%make a copy of A_prod_dissect
A0=A_prod_dissect;
Bcopy=A0;
k

%do strong RRQR for each subdomain plus the interfaces
% A_prod_dissect=full(A_prod_dissect);
last_k_columns_matrix = [];%;zeros(number_subdomains,k);
for i=1:number_subdomains
%     i=i+1
%     last_k_columns_matrix(i,:)=last_k_selected_QRCP(A0(:,A1_global_indices_matrix{i}),A0(:,A2_global_indices_matrix{i}),A0,interface_global_subdomain{i},k);
%         last_k_columns_matrix(i,:)=last_k_selected_QRCP1(A0(:,A1_global_indices_matrix{i}),A0(:,A2_global_indices_matrix{i}),A0,k);
        last_k_columns_matrix=[last_k_columns_matrix; last_k_selected_QRCP2(A0(:,A1_global_indices_matrix{i}),A0,k,k_df)'];
        last_k_columns_matrix=[last_k_columns_matrix ;last_k_selected_QRCP2(A0(:,A2_global_indices_matrix{i}),A0,k,k_df)'];
    %remove selected columns from previous iteration   
%     [A1_global_indices_matrix,A2_global_indices_matrix,A_prod_dissect,interface_global_subdomain] = RemoveSelectedcolumns(A1_global_indices_matrix,A2_global_indices_matrix,A_prod_dissect,interface_global_subdomain,last_k_columns_matrix,i,number_subdomains);
end

% figure
% spy(A_prod_dissect(:,A1_global_indices_matrix{i}))
% figure
% spy(A_prod_dissect(:,A2_global_indices_matrix{i}))

% last_k_columns_matrix

if strcmp(type,'BinaryTree')
%% Binary tree
last_k_columns_matrix_i=cell(number_subdomains,log2(npes));
for i=1:number_subdomains
    last_k_columns_matrix_i{i,1}=last_k_columns_matrix(i,:);
end

kk=1;
% ii=1;
for j=2:log2(npes)
%     j=j+1
    for i=2*kk:2*kk:number_subdomains  
%         i=i+2
        if size(last_k_columns_matrix_i{i,j},2)==0 && size(last_k_columns_matrix_i{i-1,j},2)==0 && j~=log2(npes)
          last_k_columns_matrix_i{i,j} = last_k_selected_QRCP1(A0(:,last_k_columns_matrix_i{i,j-1}),A0(:,last_k_columns_matrix_i{i-kk,j-1}),A0,k)';        
%          last_k_columns_matrix_i =[last_k_columns_matrix_i;last_k_columns_matrix_i_temp']
%         kk=kk+2;
        elseif j==log2(npes) && i== number_subdomains
            last_k_columns_matrix_i{end,end} = last_k_selected_QRCP1(A0(:,last_k_columns_matrix_i{i,j-1}),A0(:,last_k_columns_matrix_i{i-number_subdomains/2,j-1}),A0,k)';
        end
    end
    kk=kk*2;
%     ii=ii*2;
%     last_k_columns_matrix_i =[last_k_columns_matrix_i;last_k_columns_matrix_i_temp]
end
last_k_columns_matrix_i;

rejected_columns = setdiff(1:size(A,2),last_k_columns_matrix_i{end,end});
A_final = [A0(:,rejected_columns) A0(:,last_k_columns_matrix_i{end,end})];
elseif strcmp(type,'kN')
%% QRCP for k*number_subdomains columns
% %remove duplicated columns
last_k_columns_matrix_i=last_k_columns_matrix(1,:);
last_k_columns_matrix_temp = cell(number_subdomains,1);
last_k_columns_matrix_temp{1} = last_k_columns_matrix(1,:);
for i=2:2*number_subdomains
%     i=i+1
    last_k_columns_matrix_i = setdiff(last_k_columns_matrix(i,:),cat(2,last_k_columns_matrix_temp{1:i-1,:}));
    last_k_columns_matrix_temp{i} = last_k_columns_matrix_i;
end


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

% %add global interfaces
for i=1:length(global_interfaces_columns)
    BB =  [BB A0(:,global_interfaces_columns_temp{i})];
end


% BB=full(BB);
% do strong RRQR for BB
% [Q,R,P_BB] = qr(full(BB),'vector');
% k = 200;
f = 2;
[Q, R, P_BB] = sRRQR(full(BB), f, 'rank', k);


% fprintf('Step 2 selected columns ')
[~,last_k_columns_matrix_i]=ismember(BB(:,P_BB(end-k_df+1:end)).',A0.','rows');

% last_k_columns_matrix_i'

rejected_columns = setdiff(1:size(Bcopy,2),last_k_columns_matrix_i);
A_final = [Bcopy(:,rejected_columns) Bcopy(:,last_k_columns_matrix_i)];
end
P_kND = [rejected_columns last_k_columns_matrix_i'];



% p(end-k+1:end)
% p(:,last_k_columns_matrix_i)
% bb=unique([last_k_columns_matrix_i' ,rejected_columns]);
% size(bb)
% size([last_k_columns_matrix_i' ,rejected_columns])
%        
%% Check if there are overlapped columns in rejected columns
% count=0;
% for i=1:length(last_k_columns_matrix_i) %11780 %size(A)=11790=size(A_final)
%     for j=1:length(rejected_columns)
%         if last_k_columns_matrix_i(i)==rejected_columns(j)
%             count=count+1;
%         end
%     end
% end
% count
% return
%%
% figure
% spy(A_prod_dissect(A1_global_indices_matrix{1},A1_global_indices_matrix{1}))
% figure
% spy(A_prod_dissect(:,A1_global_indices_matrix{1}))
% figure
% spy(A_prod_dissect(A2_global_indices_matrix{1},A2_global_indices_matrix{1}))



%%
%----------------------------------------------------
% [Q_final,R_final] = qr(A_final,0);
%----------------------------------------------------


% [Q_final,R_final] = qr(full(A_final),0);
% [Q_final1,R_final1] = qr(full(A_final),0);
% 
% [~,S,~] = svds(A_final);
% [~,S1,~] = svd(full(A_final));
% 
% [~,S2,~] = svds(A);
% [~,S3,~] = svd(full(A));

% return
% norm(A_final,'fro')%7.3930e+04
% norm(full(A_final),2)%5.4178e+04
% A_prod_dissect=full(A_prod_dissect);
% R_final = full(R_final);
% A=full(A);
% figure
% spy(R_final)
% spy(R_final1)
% spy(R_final(end-k+1:end,end-k+1:end))
% spy(A_final)
% A_prod_dissect=full(A_prod_dissect);
% norm(A_prod) % 0.06640760279223848
% norm(R_final(end-k+1:end,end-k+1:end),2) % 2.8919558392234524e-9
% svd(A_final)
% norm(A_prod_dissect,'fro')
% norm(R_final(end-k+1:end,end-k+1:end),'fro')
% k approximates smallest singular values by QRCP

% svd for R_final
% R_final= full(R_final);
% A_prod_dissect = full(A0);
% [U_k,Sigma_k,V_k] = svd(full(R_final(:,end-k+1:end)));




%--------------------------------------------------------------------
% [U_k,Sigma_k,V_k] = svd(full(R_final));
%--------------------------------------------------------------------





% k_smallest_sv_approx = diag(Sigma_k);

% k_smallest_sv_approx_matrix(iter,:) = k_smallest_sv_approx;
% npes=2*npes;
% end


% 
% %svd of A_prod_dissect
% A = full(A);
% [U,S,V] = svd(A);
% k_smallest_sv_svd = abs(diag(S(end-k+1:end,end-k+1:end)));
% 
% 
% % sigma(end-k)
% 
% norm(k_smallest_sv_svd,'fro')
% norm(k_smallest_sv_approx,'fro')
% 
% 
% % figure
% % plot(1:k,k_smallest_sv_svd,'b',1:k,k_smallest_sv_approx,'r')
% % % plot([k_smallest_sv_svd R_final(end-k+1:end,end-k+1:end)])
% % legend('SV by SVD','Approximate SV by QRCP on subdomains')
% 
% 
% % % %QRCP for A
% A=full(A);
% [QA,RA,PA] = qr(A,0);
% k_selected_columns_A = PA(end-k+1:end);
% 
% % [~,selected_columns_A]=ismember(A(:,PA(end-k+1:end)).',A.','rows');
% 
% rejected_columns_A = setdiff(1:size(A,2),k_selected_columns_A);
% 
% A_permuted = [A(:,rejected_columns_A) A(:,k_selected_columns_A)];
% %           A(:,[rejected_columns_A k_selected_columns_A])
% 
% 
% %QR for A_permuted
% % A_permuted=full(A_permuted);
% [Q_A_permuted,R_A_permuted] = qr(A_permuted,0);
% % k_smallest_sv_QRCP_A_permuted = abs(diag(R_A_permuted(end-k+1:end,end-k+1:end)));
% 
% [~,S_A,~] = svd(R_A_permuted(end-k+1:end,end-k+1:end));
% k_smallest_sv_QRCP_A_permuted = diag(S_A);
% 
% cmap = hsv(number_iter); % color map
% figure
% semilogy(1:k,k_smallest_sv_svd,'-db',1:k,k_smallest_sv_QRCP_A_permuted,'-^k')
% % legend('SV by SVD(A)','Approximate SV by QRCP on A')
% legend_string={'SV by SVD(A)','Approximate SV by QRCP on A'};
% % plot_vect=[];
% for i=1:number_iter
%     hold on
%   plot(1:k,k_smallest_sv_approx_matrix(i,:),'Color',cmap(i,:),'Marker','s') 
% %   legend_string=[legend_string ['Approximate SV by QRCP Nested Dissection, npes=',num2str(npes_vect(i))]]
% legend_string{i+2} = ['Approximate SV by QRCP Nested Dissection, npes=',num2str(npes_vect(i))];
% end
% legend('SV by SVD(A)','Approximate SV by QRCP on A',legend_string);
% % legend('SV by SVD(A)','Approximate SV by QRCP on A',legend_string(i))
% title(['Smallest singular values comparison for ',filename])
%  xlabel('k')
%  
% 
%  
% figure
% for i=1:number_iter
%     hold on
%   plot(1:k,k_smallest_sv_svd'./k_smallest_sv_approx_matrix(i,:),'Color',cmap(i,:),'Marker','o')
% %   legend_string=[legend_string ['Approximate SV by QRCP Nested Dissection, npes=',num2str(npes_vect(i))]]
% legend_string1{i} = ['npes=',num2str(npes_vect(i))];
% end
% legend(legend_string1);
% % legend('SV by SVD(A)','Approximate SV by QRCP on A',legend_string(i))
% title(['Ratio between svd(A) and approximate SV by QRCP Nested Dissection for ',filename])
%  xlabel('k')
%  
%  figure
% for i=1:number_iter
%     hold on
%   plot(1:k,k_smallest_sv_QRCP_A_permuted'./k_smallest_sv_approx_matrix(i,:),'Color',cmap(i,:),'Marker','x')
% %   legend_string=[legend_string ['Approximate SV by QRCP Nested Dissection, npes=',num2str(npes_vect(i))]]
% legend_string2{i} = ['npes=',num2str(npes_vect(i))];
% end
% legend(legend_string2);
% % legend('SV by SVD(A)','Approximate SV by QRCP on A',legend_string(i))
% title(['Ratio between aprroximate SV by QRCP(A) and approximate SV by QRCP Nested Dissection for ',filename])
%  xlabel('k')
%  
%  ylim([0 1.2])
%  
% %  yyaxis right
% 
% %     figure
% % %      plot(1:k,k_smallest_sv_svd,'-db',1:k,k_smallest_sv_QRCP_A_permuted,'-^k',1:k,k_smallest_sv_approx,'-sr')
% % 
% % %     yyaxis right
% % plot(1:k,k_smallest_sv_svd./k_smallest_sv_approx,'-ob',1:k,k_smallest_sv_QRCP_A_permuted./k_smallest_sv_approx,'-ok')
% % legend('SV by SVD(A)','Approximate SV by QRCP on A','Approximate SV by QRCP Nested Dissection','approx \sigma_i(A)/\sigma_i(A_k)','\sigma_i(A)/\sigma_i(A_k)')
% % % ylim([0 1])
% % 
% %     title(['Smallest singular values comparison, npes=',num2str(npes)])
% %     xlabel('k')
%%
%%%--------------------------------------------------
% A_final*P_C = Q_final*R_final=Q_final * [R11  R12]
%                                         [     R22]
% [Q_final,R_final,P_C] = qr(full(A_final),'matrix');   %   [Q,R,E] = QR(A) produces unitary Q, upper triangular R and a
%                                                 %   permutation matrix E so that A*E = Q*R. The column permutation E is
%                                                 %   chosen to reduce fill-in in R.
% %---------------------------------------------------- 
A_ND = A_prod_dissect;
A_permutted = A_final;

% svdA=svd(full(A));
% svdA_kND = svd(full(A_permutted));

% normAPC_QRnormest=normest(A_final*P_C-Q_final*R_final)
% Q1=Q_final(1:n,1:n-k);
% Q2=Q_final(1:n,end-k+1:end);
% %%Approximate for the right null space of A
% R11=R_final(1:n-k,1:n-k);
% R12=R_final(1:n-k,n-k+1:end);
% R22=R_final(end-k+1:end,end-k+1:end);
% V2_tilde = P_C*[-(R11^-1)*R12;eye(k)];
% %-------------------------------------------------------



% normR22=normest(R22) %normR22 =  3.0113e-07

% P_C_temp = eye(n,n);
% P_C_mat = P_C_temp(:,P_C);

% P_k=P_k(end-Vcols+1:end,:);
% size([-(R11^-1)*R12;eye(k)])
% size(P_k) .% (1821,5)
        
% mm=1

% normest(P_k)
% normest(A*P_k)
% Q_final*[zeros(size(R12));full(R22)];
% normest(Q_final)
% size(R11)
% size(R12)
% size(R_final)

% %% generating H
%             x_0 = zeros(n, 1);
%             x_1 = x_0;
%             r_0 = b - A_final * x_0;
%             normr = norm(r_0);
%             V(:, 1) = r_0/normr;
%             e_1 = zeros(maxit + 1, 1);
%             e_1(1) = 1;
%             c = normr * e_1;
%             for j = 1 : maxit  %perform maxit steps of gmres, solving min||c-H_my||_2 for y 
% %                 w = A_final * blcLS(Pr, begin_in, end_in, V(:, j)); % and generating V_{m+1} and H_m
%                 w = A_final * V(:, j); % and generating V_{m+1} and H_m
%                 temp = V(:, 1 : j)' * w;
%                 w = w - V(:, 1 : j) * temp;
%                 temp1 = V(:, 1 : j)' * w;
%                 w = w - V(:, 1 : j) * temp1;
%                 H(1 : j, j) = temp + temp1;
%                 H(j + 1, j) = norm(w);
%                 V(:, j + 1) = w / H(j + 1, j);
% %                 [QH, RH] = qr(H(1 : j + 1, 1 : j));
% %                 e_1 = zeros(j + 1, 1);
% %                 e_1(1) = 1;
% %                 c = normr * e_1;
% %                 res = QH' * c;
% %                 fprintf('Res(%d) = %e\n', j, abs(res(end)));
% %                 resvec{matInd, methodInd} = [resvec{matInd, methodInd}, abs(res(end))];
%             end
% %             k = maxdef;
% %             y = RH(1 : end - 1, :) \ res(1 : end - 1);
% %             x_1 = x_0 + V(:, 1 : maxit) * y;
% %             r_1 = V * (c - H * y);
% %% approximate V_k
%    [U_k, Sigma_k, P_k] = svd(H,0);