function [W] = Gram_Schmidt(U)
% U=[1 1 0 0 2; 1 0 1 0 6;-1 0 0 7 1;1 2 2 1 5;2 1 3 4 6];
n=size(U,1);
R=Gauss(U);
if(rank(R)<n)
%     disp('Cac vector khong doc lap tuyen tinh')
    disp('The vectors are not linearly independent')
else
%     disp('Cac vector doc lap tuyen tinh')
    V=zeros(size(U));
     V(1,:)=U(1,:);
    for i=2:n
        s=0;
        for j=1:i-1
            s = s + (TichVH(U(i,:),V(j,:))/TichVH(V(j,:),V(j,:)))*V(j,:);
        end
        V(i,:) = U(i,:) - s;
    end
    W=zeros(size(U));
    for i=1:n
        W(i,:)= V(i,:)/dolon(V(i,:));
    end
%     disp('Ma tran truc chuan cua A la: ');
%     disp(W);
end

