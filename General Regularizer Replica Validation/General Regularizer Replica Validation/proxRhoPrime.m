function [y] = proxRhoPrime(x,c)
%Holds the form of the prox function (computed analytically)

y = 1./(1+c)*ones(size(x));


end

