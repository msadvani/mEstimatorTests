function [y] = IQfunc(fNoise,x,L0,N0,dblErr)
%returns q*I(q), useful in solving for q0

y = zeros(size(x));
%preprocessing
%y(x<0)=zeros(size(x(x<0)));

y(x>=0) = x(x>=0).*calculateInfoVec(fNoise,x(x>=0),L0,N0,dblErr);



end

