function uC = restriction(uF,nx_coarse,nx_fine)

% nt_F = size(uF,1);
% uC = zeros(nt_F,nx_coarse+1);


    uC= uF(1,1:nx_fine/nx_coarse:end);
