function [ out] = gridMin(fun,a,b,tol)
%Attempts to find the smalles value of a function between a and b using a
%grid search

N=100;




while (abs(a-b)>tol)
    xSet = linspace(a,b,N);
    funSet = zeros(size(xSet));
    for cnt = 1:N
        x = xSet(cnt);
        funSet(cnt) = fun(x);

    end

    [val,arg]=min(funSet);

    a = xSet(arg-1)
    b = xSet(arg+1)
end



out = mean([a,b]);


end

