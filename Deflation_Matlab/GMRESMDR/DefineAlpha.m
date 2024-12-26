%%%Author: Sophie Moufawad%%%
%%%Created: 24 November 2011%%% 

function [p1, edgecut, s, beginIn, endIn] = DefineAlpha (A, part)
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
M = A+A';
M = M - diag(diag(M));
cc = length(A);
if (cc==part)
 for i = 1:part
    beginIn(i) = i;
 end
endIn = beginIn;
p1 = beginIn;
edgecut = 0;
s=0;
else
    [parti,edgecut]=metismex('PartGraphKway',M,part, 0);
%     [p,ip,sizes] = metismex('NodeNDP',A_prod,npes); % 
    xx = 1;
    for i = 1:part
      bb = find(parti == i-1);
      p1(xx:xx+length(bb)-1) = bb;
      beginIn(i) = xx;
      s(i) = length(bb);
      xx = xx+length(bb);
      endIn(i)= xx-1;
    end
     
end 

% % s
% % beginIn
% % endIn

end

