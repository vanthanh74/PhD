function uF_seq =  FineSolver(A_F,m,nt_interval,nx_fine,F)

% clear 
% clc

% nx_fine = 10
u0 = F(1:nx_fine-1,1);
uF_seq=[u0];
i_temp = 1;
for i = 2:nt_interval+1
%     i = i + 1
%     u_temp =[];
    for j = 1:m
%         j = j+1
        u_temp = (A_F\u0) + F(i_temp*(nx_fine-1)+1:(i_temp+1)*(nx_fine-1));
        u0 = u_temp;
        i_temp =  i_temp +1;        
    end
        uF_seq =[uF_seq;u_temp];
              
end