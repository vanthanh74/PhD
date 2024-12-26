function y = blcLS(Pr, begin_in, end_in, x)
%block linear solver  Pr^-1 * x  with k-way partition
n = length(begin_in);
for i = 1 : n
         y(begin_in(i) : end_in(i), :) = Pr{i}\ x(begin_in(i) : end_in(i), :);
end
end




% for i = 1 : n
% %     if size(Pr{i},1)~=0
%          y(begin_in(i) : end_in(i), :) = Pr{i}\ x(begin_in(i) : end_in(i), :);%pinv(full(Pr{i}))* x(begin_in(i) : end_in(i), :);%
% %          y(begin_in(i) : end_in(i), :) = pinv(full(Pr{i}))* x(begin_in(i) : end_in(i), :);%Pr{i}\ x(begin_in(i) : end_in(i), :);%
% % 
% %     end
% end
% %y=diag(Pr)\x;
% end