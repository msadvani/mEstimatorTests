function [out] = gridMinSearchNonVec(func, a, b, N, tol)
% same as GRIDMINSEARCH except this runs a bit slower but never needs to
% access the function func as a vector.

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
       out = gridMinSearchNonVec(func, a,a+2*(b-a),2*N,tol);
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