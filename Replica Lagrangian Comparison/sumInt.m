function [ out ] = sumInt(func, a, b, N, tol)
%Approximates an integral by summing a series, tol is the error allowed on
%doubling of grid spacing.

maxIter = 10^6; %Max number of resizings allowed before error

xSet1 = linspace(a,b,N);
I1 = sum(func(xSet1))*(b-a)/N;

xSet2 = linspace(a,b,2*N);
I2 = sum(func(xSet2))*(b-a)/(2*N);

err = abs(I2-I1);

cnt=1;

while(err>tol)
    if (cnt>maxIter)
       error('max Number of Iterations exceeded - Perhaps the Integral is Not CONVERGING...') 
    end
    N=3*N;
    xSet1 = linspace(a,b,N);
    I1 = sum(func(xSet1))*(b-a)/N;

    xSet2 = linspace(a,b,2*N);
    I2 = sum(func(xSet2))*(b-a)/(2*N);

    err = abs(I2-I1);
    cnt=cnt+1;
    
end


if (I2==0)
    disp('Warning- It is possible that the Measure has not been picked up by the grid size, try LARGER N')
    
end

out = I2;



end

