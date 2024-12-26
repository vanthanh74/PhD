%%%Author: Sophie Moufawad%%%
%%%Created: 24 November 2011%%% 

function [p1, beginIn, endIn,sizes] = DefineAlphaND (A, npes)
%DefineAlpha performs kway partitioning using metis and defines the parts for each core

%This function takes as input:

% 1.   A is a sparse matrix (A=sparse([],[],[],49,49)) for which we want to
%        solve the system Ax = b and want to find a preconditioner based on
%        ILU0 factorization. 
% 2.   part is the number of blocks that we want to have, which should be
%           related to the number of cores to be used in parallel.

%This function outputs:

% 1.    p1 the permutation vector of metismex(kway)
% 2.    edgecut the number of edges that are split between different
%       domains(output of metismex)
% 3.    s contains the sizes of each of the subdomains
% 4.    beginIn the start indices of the subdomains
% 5.    endIn the end indices of the subdomains
% M = A+A';
% M = M - diag(diag(M));
M=A'*A;
% cc = length(A);
% if (cc==part)
%  for i = 1:part
%     beginIn(i) = i;
%  end
% endIn = beginIn;
% p1 = beginIn;
% edgecut = 0;
% s=0;
% else
%     [parti,edgecut]=metismex('PartGraphKway',M,part, 0);
    [p,ip,sizes] = metismex('NodeNDP',M,npes); % 
%     xx = 1;
%     for i = 1:part
%       bb = find(parti == i-1);
%       p1(xx:xx+length(bb)-1) = bb;
%       beginIn(i) = xx;
%       s(i) = length(bb);
%       xx = xx+length(bb);
%       endIn(i)= xx-1;
%     end
%      
% end 

% % s
% % beginIn
% % endIn
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
            subdomains_sizes_matrix(npes/size(separators{i-1},2):npes/size(separators{i-1},2):size(subdomains_sizes_matrix,1),rank(subdomains_sizes_matrix)+1)= separators{i}';                      
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

subdomains_sizes_matrix;






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
subdomain_indices;
countNull=0;
for i=1:size(subdomain_indices,1)
%     i=i+1
    for j=1:size(subdomain_indices,2)
%         j=j+1
        if size(subdomain_indices{i,j},1)~=0 && size(subdomain_indices{i,j},2)~=0
            beginIn((i-1)*i+j+countNull) =  subdomain_indices{i,j}(1);
            endIn((i-1)*i+j+countNull)   =  subdomain_indices{i,j}(end);
        else
            countNull=countNull+1;
        end
    end
end
p1=p;
end

