%Did a comparison with CompareBayesOptL2 - this is the same basic code, but
%for a single kappa

clear all;
close all;

kappa = 1.325;


%fNoise = @(x) (2*pi)^(-1/2)*exp(-x.^2/2);
%fNoise = @(x) probDist(x);
fNoise = @(x) (1/2)*exp(-abs(x))
%fNoise = @(x) exp(x)./((1+exp(x)).^2);

g = @(x) (1/2)*exp(-abs(x));


% 
% %% Theory for Optimal error
% 
% L0 = 10; %initial interval length
% N0=10^3; %initial num point in integral approx
% dblErr = .0001; %permissible error in I
% 
% 
% I = @(x) calculateInfoMat(fNoise,x,L0,N0, dblErr);
% J = @(x) calculateInfoMat(g,x, L0, N0, dblErr); 
% 
% 
% 
% %plotting the function
% %a = linspace(0,50,200);
% %plot(a,a - a.^2.*J(a))
% 
% %two constraints
% F1 = @(q,a) (a - a.^2.*J(a) - q);
% F2 = @(q,a) (I(q).*a - kappa);
% 
% 
% 
% % numA = 50;
% % numQ = 50;
% % aSet =repmat(linspace(0,400,numA),numQ,1);
% % qSet = repmat(linspace(0,10,numQ)',1,numA);
% 
% 
% %aMin1 = @(q)gridMinSearch(@(x)F1(q,x),0,5,10,.01);
% %aMin2 = @(q)gridMinSearch(@(x)F2(q,x),0,5,10,.01);
% 
% 
% aMin1 = @(q)findZeroBB(@(x)F1(q,x),.1,10,.001)
% aMin2 = @(q)findZeroBB(@(x)F2(q,x),.1,10,.001)
% 
% %d = @(q)(aMin1(q)-aMin2(q)).^2;
% %qOpt = gridMinSearchNonVec(d,.05,2,10,.005)
% 
% d = @(q)(aMin1(q)-aMin2(q));
% qOpt = findZeroBB(d,.05,2,.001);
% 
% a1 = aMin1(qOpt);
% a2 = aMin2(qOpt);
% 
% aOpt = .5*(a1+a2);
% 

%qOpt
%aOpt

qOpt = .01;
aOpt = 2;

gauss = @(x,sigma) (2*pi*sigma.^2).^(-1/2)*exp(-x.^2./(2*sigma.^2));
gaussPrime = @(x,sigma) -(2*pi*sigma.^6).^(-1/2).*x.*exp(-x.^2./(2*sigma.^2));


zeta = @(x)convolNumeric(x,@(y)fNoise(y), @(y)gauss(y,sqrt(qOpt)));
zetaPrime = @(x)convolNumeric(x,@(y)fNoise(y), @(y)gaussPrime(y,sqrt(qOpt)));

%zeta2 = @(x)convolNumeric(x,@(y)gauss(y,sqrt(qOpt)),@(y)fNoise(y)); %just a check
xi = @(x)convolNumeric(x,@(y)g(y), @(y)gauss(y,sqrt(aOpt)));

p2 = @(x) x.^2/2;


%r = @(x)(p2(x)+qOpt*ln(zeta(x)));

R = @(x) (p2(x) + qOpt*log(zeta(x)));  %This is the function you want to conjugate

Rprime = @(x) (x+qOpt*zetaPrime(x)./zeta(x))

%Validation of Rprime
% eps=.0001;
% RprimeApprox = @(x)(R(x+eps)-R(x))/eps
% 
% 
% 
% %plot(xSet,-log(xi(xSet)))
% plot(xSet,RprimeApprox(xSet),'ko')

numX =30;
xSet =linspace(-6,6,numX);
hold on;
%plot(xSet,R(xSet))
%plot(xSet,Rprime(xSet),'r')

FCform = @(x,y)(x.*y-R(y))
FCformPrime = @(x,y)(x-Rprime(y))
FCformPrimeNeg = @(x,y)(Rprime(y)-x)
%plot(xSet,FCform(-2,xSet))
%plot(xSet,FCformPrime(-2,xSet),'r')



%Confirmation of FCformPrime
% eps = .0001;
% FCformPrimeApprox = @(x,y)((FCform(x,y+eps)-FCform(x,y))/eps);
% plot(xSet,FCformPrimeApprox(2,xSet),'g*')



ySoln = @(x)findZeroBB_realLine(@(y)FCformPrimeNeg(x,y),.1,10,.001);
conj = @(x) FCform(x,ySoln(x));

% conjSet = zeros(size(xSet));
% for cnt = 1:length(xSet)
%     [cnt, length(xSet)]
%    conjSet(cnt) = conj(xSet(cnt));
% end
% 
% plot(xSet,conjSet,'k--')


rhoOpt = @(x)(conj(x)-p2(x))


rhoOptSet = zeros(size(xSet));
for cnt = 1:length(xSet)
    rhoOptSet(cnt) = rhoOpt(xSet(cnt));
end
 
plot(xSet,rhoOptSet,'k--')


%Now you need to optimize this by finding ySoln(x)






%plot(xSet,-log(fNoise(xSet)))



