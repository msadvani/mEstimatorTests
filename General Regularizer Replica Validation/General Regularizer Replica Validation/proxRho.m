function [y] = proxRho(x,c)
%Holds the form of the prox function (computed analytically)

y = x./(1+c);


end

