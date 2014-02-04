function [qOut,cOut] = qThyL2(k, lambda, varCoeff, varNoise)
%Exact analytic error for two l2 norms with a lambda ratio
%TO DO: add in additional parameters dependent on the variances of f and g!

b = 1-k-k.*lambda;
c = -k.*lambda;
 
d = (-b+sqrt(b.^2-4*c))/2;

 
qOut = ((varCoeff*d.^2 + varNoise*k)./((1+d).^2 - k));

cOut = d/(lambda*k) - 1;


end

