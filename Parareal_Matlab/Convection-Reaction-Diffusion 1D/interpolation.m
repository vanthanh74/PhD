function uF = interpolation(uC,nx_fine)

nt_C = size(uC,1);
uF = zeros(nt_C,nx_fine+1);

for i =1:nt_C
    uF(i,1:2:end) = uC(i,:);
    uF(i,2:2:end-1) = ( uF(i,1:2:end-1)+uF(i,3:2:end) )/2;
end