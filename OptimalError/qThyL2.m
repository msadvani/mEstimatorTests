function [out] = qThyL2(k, lambda)
%Exact analytic error for two l2 norms with a lambda ratio

b = 1-k-k.*lambda;
c = -k.*lambda;
 
d = (-b+sqrt(b.^2-4*c))/2;

 
out = ((d.^2 + k*2)./((1+d).^2 - k));



end

