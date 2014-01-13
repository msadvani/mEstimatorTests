function [y] = probDist(x)
%Defines a probability distribution of interest

y = zeros(size(x));

%Exponential
y(x>0) = exp(-x(x>0));
y(x<0) = zeros(size(x(x<0)));


end

