function [y] = proxSigma(x,c)
%Holds the form of the prox function (computed analytically)

%lambda =1 L2 norm regularizer
%y = x./(1+c);

y = zeros(size(x));
y(abs(x)<c) = zeros(size(find(abs(x)<c)));
y(x>c) = x(x>c)-c;
y(x<-c) = x(x<-c)+c;



end

