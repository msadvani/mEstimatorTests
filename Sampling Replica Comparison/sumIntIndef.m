function [ out ] = sumIntIndef(func, L, N, tol)
%Approximates an integral by summing a series, tol is the error allowed on
%doubling of grid spacing, and also resizing. The integral is between -L
%and L

maxIter = 10^6; %Max number of resizings allowed before error

xSet1 = linspace(-L,L,N);
I1 = sum(func(xSet1))*(2*L)/N;

xSet2 = linspace(-L,L,2*N);
I2 = sum(func(xSet2))*(2*L)/(2*N);

xSet3 = linspace(-2*L,2*L,2*N);
I3 = sum(func(xSet3))*(4*L)/(2*N);



err12 = abs(I2-I1); %error from doubling grid 
err13 = abs(I3-I1); %error from doubling interval

cnt=1;


while(err12>tol || err13>tol)
    if (cnt>maxIter)
       error('max Number of Iterations exceeded - Perhaps the Integral is Not CONVERGING...') 
    end
    
    
    if (err12>tol)
        N=3*N;
    elseif(err13>tol)
        L=3*L;
        N=3*N;
    end
        
    xSet1 = linspace(-L,L,N);
    I1 = sum(func(xSet1))*(2*L)/N;

    xSet2 = linspace(-L,L,2*N);
    I2 = sum(func(xSet2))*(2*L)/(2*N);

    xSet3 = linspace(-2*L,2*L,2*N);
    I3 = sum(func(xSet3))*(4*L)/(2*N);



    err12 = abs(I2-I1); %error from doubling grid 
    err13 = abs(I3-I1); %error from doubling interval

end


if (I3==0)
    disp('Warning- It is possible that the Measure has not been picked up by the grid size, try LARGER N')
    
end

out = I3;

end

