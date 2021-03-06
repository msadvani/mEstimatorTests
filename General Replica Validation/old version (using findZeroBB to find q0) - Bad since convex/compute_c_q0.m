clear all;
close all;

kappa = .8;

gauss = @(x,sigma) (2*pi*sigma.^2).^(-1/2)*exp(-x.^2./(2*sigma.^2));

sigNoise = 2;
fNoise = @(x) gauss(x,sigNoise);
g = @(x) 1/2*exp(-abs(x));

varNoise = sumIntIndef(@(x) x.^2.*fNoise(x), 2,10,.00001);
varCoeff = sumIntIndef(@(x) x.^2.*g(x), 2,10,.00001)

c0Hat = @(q0,c) (1/(2*kappa*(1+c).^2))*(q0 + varNoise);
cHat = @(q0,c) (1./(2*kappa*(1+c)));

a = @(q0) (q0 + varNoise)*kappa;


N=5000;
L=13;







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


[q0Thy,cThy] = qThyL2(kappa, 1, varCoeff, varNoise)
%Now we need to compare this with the analytic solution already worked out!

%computeConstraint1(1,1,fNoise, g,kappa,N,L,.1)
F1 = @(q0,c)computeConstraint1(q0,c,fNoise, g,kappa,N,L,.01)
F2 = @(q0,c)computeConstraint2(q0,c,fNoise, g,kappa,N,L,.01)



lambda =1;
F1Thy = @(q0,c) -q0 + (lambda^2*kappa^2*(1+c).^2*varCoeff + a(q0))./(1+lambda*kappa*(1+c)).^2
%F1Thy = @(q0,c) q0.*(1+lambda*kappa*(1+c)).^2 - (lambda^2*kappa^2*(1+c).^2*varCoeff + a(q0))

%F1Thy = @(q0,c) q0 - 1/(1+kappa*(1+c))^2*(kappa^2*(1+c)^2*varCoeff - a(q0,c)) 

F2Thy = @(q0,c) c./(1+c) - kappa./(1+lambda*kappa*(1+c))





%for plotting
q0=1.46;



numC=100;
cSet = linspace(0,10,numC);
 
F1setThy = zeros(size(cSet));
F1set = zeros(size(cSet));
 

F2setThy = zeros(size(cSet));
F2set = zeros(size(cSet));
 
for cnt =1:numC
    [cnt,numC] 
    c = cSet(cnt);
    F1set(cnt) = F1(q0,c);
    F1setThy(cnt) = F1Thy(q0,c);
    F2set(cnt) = F2(q0,c);
    F2setThy(cnt) = F2Thy(q0,c);

end
 
 hold on;

plot(cSet,F1set,'-')
plot(cSet,F1setThy,'r*')

plot(cSet,F2set,'k')
plot(cSet,F2setThy,'g*')

plot(cSet,zeros(size(cSet)),'k')




tic

c1Min = @(q) findZeroBB_cutoff(@(c)F1(q,c),.01,2,.001);
c2Min = @(q) findZeroBB_cutoff(@(c)F2(q,c),.01,2,.001);









%qOpt = findZeroBB(@(q0) (c1Min(q)-c2Min(q)),.01,2,.001);





% numQ = 20;
% 
% c1Set = zeros(1,numQ);
% c2Set = zeros(1,numQ);
% 
% qSet = linspace(0,2,numQ)
% for qCnt =1:numQ
%      [qCnt,numQ] 
%      q = qSet(qCnt);
%      c1Set(qCnt) = c1Min(q);
%      c2Set(qCnt) = c2Min(q);
% end
% 
% 
% hold on;
% 
% plot(qSet,c1Set,'-')
% plot(qSet,c2Set,'r-')
















%q0 = findZeroBB(@(q) (c1Min(q)-c2Min(q)),.01,2,.001)


cVec=[c1Min(q0),c2Min(q0)]

toc





































