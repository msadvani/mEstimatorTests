function [I] = int2D(func,L0,N0,tol)
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


