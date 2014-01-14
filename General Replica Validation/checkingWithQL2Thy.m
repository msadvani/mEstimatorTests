
%f = @(x,y) (2*pi)^(-1)*exp(-x.^2/2-y.^2/2)
%intImp2D(f,5,50,.0001)

gauss = @(x,sigma) (2*pi*sigma.^2).^(-1/2)*exp(-x.^2./(2*sigma.^2));

sigNoise = 2;
fNoise = @(x) gauss(x,sigNoise);
g = @(x) 1/2*exp(-abs(x));



kappa = .5;




varNoise = sumIntIndef(@(x) x.^2.*fNoise(x), 2,10,.00001);
varCoeff = sumIntIndef(@(x) x.^2.*g(x), 2,10,.00001)


%% Numerically finding c0hat - uneccisary if we assume rho =x^2/2
%[I, L,N]=intImp2D(@(x,z)(fNoise(x).*gaussUnit(z)),10,100,.001)
%L=10,N=1600 (3 digits of accuracy .001)
%[I, L,N]=intImp2D(@(x,z)(fNoise(x).*gaussUnit(z)),10,1600,.001)
%This will be superflous if we assume rho = x^2/2
%c0Hat = @(q0,c) 1/(2*kappa*c.^2)*intImp2D(@(x,z) fNoise(x).*gaussUnit(z).*(proxRho(sqrt(q0)*z + x,c) - sqrt(q0)*z - x).^2,10,2000,.001)




c0Hat = @(q0,c) (1/(2*kappa*(1+c).^2))*(q0 + varNoise);
cHat = @(q0,c) (1./(2*kappa*(1+c)));

a = @(q0,c) (q0 + varNoise)*kappa;





N=10000;
L=12;

xSet = linspace(-L,L,N);
gSet = g(xSet);

normSet = @(q0,c) gauss(xSet,sqrt(a(q0,c)));

xi = @(q0,c) conv(normSet(q0,c),gSet)*(2*L/N);

convMeas = linspace(-2*L,2*L,2*N-1);
dConvMeas = (4*L/(2*N-1));


%Note that as q0 becomes very very large this is not longer accurate.
%Perhaps you should require choosing N, L s.t I(q0, c) >.995
I = @(q0,c) sum(xi(q0,c))*(4*L/(2*N-1))


%eventually you should define these as sepearate functions...

F1= @(q0,c) q0 + a(q0,c) - sum(xi(q0,c).*(convMeas - proxSigma(convMeas,1./(2*cHat(q0,c)))).^2*dConvMeas) - 2*a(q0,c)*sum(xi(q0,c).*proxSigmaPrime(convMeas,1./(2*cHat(q0,c)))*dConvMeas);


%This should be correct even if convMeas not right...
F2 = @(q0,c) -c./(1+c)+kappa*sum(xi(q0,c).*proxSigmaPrime(convMeas,1/(2*cHat(q0,c))))*dConvMeas;


[q0Thy,cThy] = qThyL2(kappa, 1, varCoeff, varNoise)
%Now we need to compare this with the analytic solution already worked out!









%write a convolution version...to make computation faster

%%F2 = @(q0,c) (1/(1+c)-1+kappa*intImp2D(@(x,y) proxSigma(x+sqrt(a)*y,1/(2*cHat(q0,c))),10,2000));


























