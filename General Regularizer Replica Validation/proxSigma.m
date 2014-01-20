function [y] = proxSigma(x,c,form,lambda)
%Holds the form of the prox function (computed analytically)

c=c*lambda;
if (form ==1)
    y = zeros(size(x));
    y(abs(x)<c) = zeros(size(find(abs(x)<c)));
    y(x>c) = x(x>c)-c;
    y(x<-c) = x(x<-c)+c;

elseif (form==2)  %L2 norm regularizer
    y = x./(1+c);
end

end

