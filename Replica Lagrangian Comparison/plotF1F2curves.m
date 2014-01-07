clear all;
close all;

kappa = 2.3;


fNoise = @(x) (2*pi)^(-1/2)*exp(-x.^2/2);
g = @(x) (1/2)*exp(-abs(x));


%% Theory for Optimal error

L0 = 10; %initial interval length
N0=10^3; %initial num point in integral approx
dblErr = .0001; %permissible error in I


I = @(x) calculateInfoMat(fNoise,x,L0,N0, dblErr);
J = @(x) calculateInfoMat(g,x, L0, N0, dblErr); 



%plotting the function
%a = linspace(0,50,200);
%plot(a,a - a.^2.*J(a))

%two constraints
F1 = @(q,a) (a - a.^2.*J(a) - q)
F2 = @(q,a) (I(q).*a - kappa);


q = 1; %Set q values

aSet = linspace(.1,10,100);


hold on;

plot(aSet,F1(q,aSet),'r');
plot(aSet,F2(q,aSet));
plot(aSet, zeros(size(aSet)),'k');

legend('F1','F2')

xlabel('a')
xlabel('F')



findZeroBB(@(x)F1(q,x),1,2,.01)
