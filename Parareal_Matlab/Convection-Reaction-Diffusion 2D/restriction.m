function uC = restriction(uF,nx_coarse,ny_coarse,nx_fine,ny_fine)
% nx_coarse = 16
% nx_fine = 81
% uF = u0'
% nt_F = size(uF,1);
% uC = zeros(1,nx_coarse);

% uC = zeros(nx_coarse,ny_coarse);

    uC = uF(1:ny_fine/ny_coarse:end,1:nx_fine/nx_coarse:end);


