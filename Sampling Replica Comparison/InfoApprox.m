function [I] = InfoApprox(fNoise,q0,L,N)
    sigma = sqrt(q0);
    gauss = @(x)((2*pi.*sigma^2)^(-.5)*exp(-x.^2./(2.*sigma^2)));
    gaussPrime = @(x)(-(2*pi.*sigma^2)^(-.5)*(x/(sigma^2)).*exp(-x.^2./(2.*sigma^2)));
    zSet = linspace(-L,L,N);
    dx = (2*L)/N;
    dxConv = (4*L)/(2*N-1); %delta for the convolution function
    zeta = conv(fNoise(zSet),gauss(zSet))*dx;
    zetaPrime = conv(fNoise(zSet),gaussPrime(zSet))*dx;
    Ifunc = (zetaPrime.^2)./zeta;
    Ifunc(isnan(Ifunc))=0;
    I = sum(Ifunc)*dxConv;
end

