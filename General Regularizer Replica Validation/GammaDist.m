function [y] = GammaDist(x)
%GAMMADISTRIB currently only outputs the simplest gamma distrib
%values...can be updated soon


y = zeros(size(x));


y(x>0) = x(x>0).*exp(-x(x>0));





end

