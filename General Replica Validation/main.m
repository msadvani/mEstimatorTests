
%f = @(x,y) (2*pi)^(-1)*exp(-x.^2/2-y.^2/2)
%intImp2D(f,5,50,.0001)


fNoise = @(x) 1/2*exp(-abs(x));
gaussUnit = @(x) (2*pi)^(-1/2)*exp(-x.^2/2);


kappa = .5;

varNoise = sumIntIndef(@(x) x.^2.*fNoise(x), 2,10,.00001);
%varCoeff = sumIntIndef(@(x) x.^2.*g(x), 2,10,.00001)


%% Numerically finding c0hat - uneccisary if we assume rho =x^2/2
%[I, L,N]=intImp2D(@(x,z)(fNoise(x).*gaussUnit(z)),10,100,.001)
%L=10,N=1600 (3 digits of accuracy .001)
%[I, L,N]=intImp2D(@(x,z)(fNoise(x).*gaussUnit(z)),10,1600,.001)
%This will be superflous if we assume rho = x^2/2
%c0Hat = @(q0,c) 1/(2*kappa*c.^2)*intImp2D(@(x,z) fNoise(x).*gaussUnit(z).*(proxRho(sqrt(q0)*z + x,c) - sqrt(q0)*z - x).^2,10,2000,.001)




c0Hat = @(q0,c) (1/(2*kappa*(1+c).^2))*(q0 + varNoise);
cHat = @(q0,c) (1./(2*kappa*(1+c)));

a = @(q0,c) (q0 + varNoise)*kappa;

%write a convolution version...to make computation faster

%%F2 = @(q0,c) (1/(1+c)-1+kappa*intImp2D(@(x,y) proxSigma(x+sqrt(a)*y,1/(2*cHat(q0,c))),10,2000));


























