function [out] = gridMinSearchNonVec(func, a, b, N, tol)
%GRIDMINSEARCH. Returns the value that minimizes the function to a given
%tolerance

    xSet = linspace(a,b,N);
    funSet = zeros(size(xSet));
    for cnt = 1:N
        funSet(cnt) = func(xSet(cnt));
    end
    
    [~, arg]=min(funSet);
    %out = xSet(arg);
    
    
    if ((b-a)/N)<tol
       out = xSet(arg);
    elseif (arg==N)
       out = gridMinSearchNonVec(func, a,a+2*(b-a),N,tol);
    else
       if(arg==1)
           a2=xSet(1);
       else
           a2= xSet(arg-1);
       end
       b2 = xSet(arg+1);
       out = gridMinSearchNonVec(func, a2, b2, N, tol);
    end
    
end