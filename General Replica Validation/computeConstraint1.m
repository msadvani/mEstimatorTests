function [ F1, I ] = computeConstraint1(q0,c, fNoise, g,kappa,N0,L0,pTol)
%Returns the value of the 1st constraint, where the integral over the
%convolved coefficient distribution is approximately 1.

gauss = @(x,sigma) (2*pi*sigma.^2).^(-1/2)*exp(-x.^2./(2*sigma.^2));

N = N0;
L = L0;

varNoise = sumIntIndef(@(x) x.^2.*fNoise(x), 2,10,.00001);

c0Hat = @(q0,c) (1/(2*kappa*(1+c).^2))*(q0 + varNoise);
cHat = @(q0,c) (1./(2*kappa*(1+c)));

a = @(q0,c) (q0 + varNoise)*kappa;


xSet = linspace(-L,L,N);
gSet = g(xSet);

normSet = @(q0,c) gauss(xSet,sqrt(a(q0,c)));

xi = @(q0,c) conv(normSet(q0,c),gSet)*(2*L/N);

convMeas = linspace(-2*L,2*L,2*N-1);
dConvMeas = (4*L/(2*N-1));


%Note that as q0 becomes very very large this is not longer accurate.
%Perhaps you should require choosing N, L s.t I(q0, c) >.995
%I = @(q0,c) sum(xi(q0,c))*(4*L/(2*N-1))
I = sum(xi(q0,c))*(4*L/(2*N-1));

if (abs(1-I)>pTol)
    error('Poor estimate of probability normalization');
end

F1= - q0 - a(q0,c) + sum(xi(q0,c).*(convMeas - proxSigma(convMeas,1./(2*cHat(q0,c)))).^2*dConvMeas) + 2*a(q0,c)*sum(xi(q0,c).*proxSigmaPrime(convMeas,1./(2*cHat(q0,c)))*dConvMeas);

% while abs(1-I)>pTol
%     L=L+L0;
%     
%     
%     xSet = linspace(-L,L,N);
%     gSet = g(xSet);
% 
%     normSet = @(q0,c) gauss(xSet,sqrt(a(q0,c)));
% 
%     xi = @(q0,c) conv(normSet(q0,c),gSet)*(2*L/N);
% 
%     convMeas = linspace(-2*L,2*L,2*N-1);
%     dConvMeas = (4*L/(2*N-1));
%     
%     I = sum(xi(q0,c))*(4*L/(2*N-1))
% end
% 
% %returns F1 only when I is close enough to 1.
% F1= q0 + a(q0,c) - sum(xi(q0,c).*(convMeas - proxSigma(convMeas,1./(2*cHat(q0,c)))).^2*dConvMeas) - 2*a(q0,c)*sum(xi(q0,c).*proxSigmaPrime(convMeas,1./(2*cHat(q0,c)))*dConvMeas);


end

