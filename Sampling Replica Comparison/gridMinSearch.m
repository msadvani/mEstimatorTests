function [out] = gridMinSearch(func, a, b, N, tol)
%GRIDMINSEARCH. Returns the value that minimizes the function to a given
%tolerance

    xSet = linspace(a,b,N);
    [val, arg]=min(func(xSet));
    
    if ((b-a)/N)<tol
       out = xSet(arg);
    elseif (arg==N)
       out = gridMinSearch(func, a,a+2*(b-a),N,tol);
    else
       a2= xSet(arg-1);
       b2 = xSet(arg+1);
       out = gridMinSearch(func, a2, b2, N, tol);
    end
    
 
end