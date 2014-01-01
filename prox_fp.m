function [y] = prox_fp(x,c)
%Derivative of a prox of a function of interest with a free parameter c

    y = zeros(size(x));
    y(x>c)=ones(size(x(x>c)));
    y(x<-c)=ones(size(x(x<-c)));
    y(x>=-c & x<=c) = zeros(size(x(x>=-c & x<=c))); 

end

