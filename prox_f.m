function [y] = prox_f(x,c)
%Prox of a function of interest with a free parameter c

    y = zeros(size(x));
    y(x>c)=x(x>c) - c;
    y(x<-c)=x(x<-c) + c;
    y(x>=-c & x<=c) = zeros(size(x(x>=-c & x<=c))); 

end

