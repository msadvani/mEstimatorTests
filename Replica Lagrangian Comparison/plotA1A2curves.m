clear all;
close all;

kappa = 50;


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
F1 = @(q,a) (a - a.^2.*J(a) - q);
F2 = @(q,a) (I(q).*a - kappa);


% numA = 50;
% numQ = 50;
% aSet =repmat(linspace(0,400,numA),numQ,1);
% qSet = repmat(linspace(0,10,numQ)',1,numA);


%aMin1 = @(q)gridMinSearch(@(x)F1(q,x),0,5,10,.01);
%aMin2 = @(q)gridMinSearch(@(x)F2(q,x),0,5,10,.01);

aMin1 = @(q)findZeroBB(@(x)F1(q,x),.1,10,.01)
aMin2 = @(q)findZeroBB(@(x)F2(q,x),.1,10,.01)



%% Now plot aMin1 and aMin2

numQ = 10;

qSet = linspace(.05,4,numQ);

aMin1Set = zeros(size(qSet));
aMin2Set = zeros(size(qSet));


for qcnt = 1:numQ
    [qcnt,numQ]
    qVal = qSet(qcnt);
    aMin1Set(qcnt) = aMin1(qVal);
    aMin2Set(qcnt) = aMin2(qVal); 
end


hold on
plot(qSet, aMin1Set,'r')
plot(qSet, aMin2Set)

legend('a1','a2')

xlabel('q')
ylabel('a')





