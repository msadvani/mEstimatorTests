clear all;
close all;

kappa = .5;

gauss = @(x,sigma) (2*pi*sigma.^2).^(-1/2)*exp(-x.^2./(2*sigma.^2));

sigNoise = 1;
fNoise = @(x) gauss(x,sigNoise);
g = @(x) 1/2*exp(-abs(x));

varNoise = sumIntIndef(@(x) x.^2.*fNoise(x), 2,10,.00001);
varCoeff = sumIntIndef(@(x) x.^2.*g(x), 2,10,.00001)

c0Hat = @(q0,c) (1/(2*kappa*(1+c).^2))*(q0 + varNoise);
cHat = @(q0,c) (1./(2*kappa*(1+c)));

a = @(q0) (q0 + varNoise)*kappa;


N=10000;
L=15;







% xSet = linspace(-L,L,N);
% gSet = g(xSet);
% 
% normSet = @(q0,c) gauss(xSet,sqrt(a(q0,c)));
% 
% xi = @(q0,c) conv(normSet(q0,c),gSet)*(2*L/N);
% 
% convMeas = linspace(-2*L,2*L,2*N-1);
% dConvMeas = (4*L/(2*N-1));
% 
% 
% %Note that as q0 becomes very very large this is not longer accurate.
% %Perhaps you should require choosing N, L s.t I(q0, c) >.995
% I = @(q0,c) sum(xi(q0,c))*(4*L/(2*N-1));
% 
% 
% %eventually you should define these as sepearate functions...
% 
% F1= @(q0,c) q0 + a(q0,c) - sum(xi(q0,c).*(convMeas - proxSigma(convMeas,1./(2*cHat(q0,c)))).^2*dConvMeas) - 2*a(q0,c)*sum(xi(q0,c).*proxSigmaPrime(convMeas,1./(2*cHat(q0,c)))*dConvMeas);
% 
% 
% %This should be correct even if convMeas not right...



F2 = @(q0,c) -c./(1+c)+kappa*sum(xi(q0,c).*proxSigmaPrime(convMeas,1/(2*cHat(q0,c))))*dConvMeas;


[q0Thy,cThy] = qThyL2(kappa, 1, varCoeff, varNoise);
%Now we need to compare this with the analytic solution already worked out!

%computeConstraint1(1,1,fNoise, g,kappa,N,L,.1)
F1 = @(q0,c)computeConstraint1(q0,c,fNoise, g,kappa,N,L,.01)
F2 = @(q0,c)computeConstraint2(q0,c,fNoise, g,kappa,N,L,.01)



lambda =1;
F1Thy = @(q0,c) -q0 + (lambda^2*kappa^2*(1+c).^2*varCoeff + a(q0))./(1+lambda*kappa*(1+c)).^2
%F1Thy = @(q0,c) q0.*(1+lambda*kappa*(1+c)).^2 - (lambda^2*kappa^2*(1+c).^2*varCoeff + a(q0))

%F1Thy = @(q0,c) q0 - 1/(1+kappa*(1+c))^2*(kappa^2*(1+c)^2*varCoeff - a(q0,c)) 



F2Thy = @(q0,c) c./(1+c) - kappa./(1+lambda*kappa*(1+c))







































