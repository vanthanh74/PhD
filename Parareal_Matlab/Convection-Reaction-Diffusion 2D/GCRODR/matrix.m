% construct the coefficient matrix
clear all
clc
load B.txt;
b=B(:,1)+B(:,2)*i;
load M_MTRIX.txt;
M1=reshape(M_MTRIX(:,3),1803,1803);
M2=reshape(M_MTRIX(:,4),1803,1803);
M=M1+M2*i;
M=M';
load N_MTRIX.txt;
N1=reshape(N_MTRIX(:,3),1803,1803);
N2=reshape(N_MTRIX(:,4),1803,1803);
N=N1+N2*i;
N=N';

