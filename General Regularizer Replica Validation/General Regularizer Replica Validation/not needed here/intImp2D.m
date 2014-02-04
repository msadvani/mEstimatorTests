function [ out, Lout,Nout] = intImp2D(func,L0,N0,tol)
%Approximates an improper integral from -inf to inf of a 2D function, given
%L0, N0 initial conditions and tolerance under doubling of integral
%sampling / size

N=N0;
L=L0;

maxIter = 10^6;  %End the function if it takes too long
numIter = 0;

xSet = repmat(linspace(-L,L,N),N,1);
ySet = repmat(linspace(-L,L,N)',1,N);

I = sum(sum(func(xSet,ySet)))*4*(L/N)^2;


xSet = repmat(linspace(-L,L,2*N),2*N,1);
ySet = repmat(linspace(-L,L,2*N)',1,2*N);

I2N = sum(sum(func(xSet,ySet)))*(L/N)^2;



xSet = repmat(linspace(-2*L,2*L,2*N),2*N,1);
ySet = repmat(linspace(-2*L,2*L,2*N)',1,2*N);

I2L = sum(sum(func(xSet,ySet)))*4*(L/N)^2;


err2N = abs(I2N-I);
err2L = abs(I2L-I);

while(err2N>tol || err2L>tol)
    numIter = numIter+1
    if(numIter>maxIter)
        error('Too many iterations to convergence')
        
    end
    if(err2N>tol)
        N = 2*N;
    elseif(err2L>tol)
        L = 2*L;
        N = 2*N;
    end
    
    
    xSet = repmat(linspace(-L,L,N),N,1);
    ySet = repmat(linspace(-L,L,N)',1,N);

    I = sum(sum(func(xSet,ySet)))*4*(L/N)^2


    xSet = repmat(linspace(-L,L,2*N),2*N,1);
    ySet = repmat(linspace(-L,L,2*N)',1,2*N);

    I2N = sum(sum(func(xSet,ySet)))*(L/N)^2


    xSet = repmat(linspace(-2*L,2*L,2*N),2*N,1);
    ySet = repmat(linspace(-2*L,2*L,2*N)',1,2*N);

    I2L = sum(sum(func(xSet,ySet)))*4*(L/N)^2

    
    err2N = abs(I2N-I);
    err2L = abs(I2L-I);
end



out = I2L;

Lout = L;
Nout = N;
