function [X]=dolon(a)
X=0;
for i=1:length(a)
    X=X + a(i)^2;
end
X=sqrt(X);