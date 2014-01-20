function [y] = proxSigmaPrime(x,c,form, lambda)
%Holds the form of the prox function (computed analytically)
c = lambda*c;

if (form==1)
    y = zeros(size(x));
    y(abs(x)>c) = ones(size(find(abs(x)>c)));
elseif (form ==2)
    y = 1./(1+c)*ones(size(x));
end

end

