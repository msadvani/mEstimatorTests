clear all;
close all;


gauss = @(x,sigma) (2*pi*sigma.^2).^(-1/2)*exp(-x.^2./(2*sigma.^2));

sigNoise = 2;
fNoise = @(x) gauss(x,sigNoise);
g = @(x) 1/2*exp(-abs(x));

computeError(.5,fNoise, g)


























