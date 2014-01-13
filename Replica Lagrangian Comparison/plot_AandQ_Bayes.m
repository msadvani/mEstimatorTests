%This code compares the optimal L2 based loss function error with the
%predicted optimal erro based on my research. qOpt should always be lower,
%apart from numerical inaccuracy.

clear all;
close all;
 
numKappa = 20;
kappaSet = linspace(.1,5,numKappa);


fNoise = @(x) (2*pi)^(-1/2)*exp(-x.^2/2); %noise
%fNoise = @(x) (1/2)*exp(-abs(x-10)); %noise
%fNoise = @(x) exp(x)./((1+exp(x)).^2);


%fNoise = @(x) probDist(x);



%g = @(x) (2*pi)^(-1/2)*exp(-x.^2/2); %coefficient distrib.
g = @(x) (1/2)*exp(-abs(x)); %coefficient distrib.

%varNoise = sumIntIndef(@(x) x.^2.*fNoise(x), 2,10,.00001)
%varCoeff = sumIntIndef(@(x) x.^2.*g(x), 2,10,.00001)




qOpt = zeros(1, numKappa);
aOpt = zeros(2, numKappa);
tic
for kcnt = 1:numKappa
    [kcnt,numKappa]

    kappa = kappaSet(kcnt);

    %% Theory for Optimal error

    L0 = 10; %initial interval length
    N0=10^3; %initial num point in integral approx
    dblErr = .0001; %permissible error in I


    I = @(x) calculateInfoMat(fNoise,x,L0,N0, dblErr);
    J = @(x) calculateInfoMat(g,x, L0, N0, dblErr);
    
    %two constraints
    F1 = @(q,a) (a - a.^2.*J(a) - q);
    F2 = @(q,a) (I(q).*a - kappa);


    aMin1 = @(q)findZeroBB(@(x)F1(q,x),.05,3,.001);
    aMin2 = @(q)findZeroBB(@(x)F2(q,x),.05,3,.001);


    d = @(q)(aMin1(q)-aMin2(q));

    qOpt(kcnt) = findZeroBB(d,.05,2,.001);
    
    aOpt(:,kcnt) = [aMin1(qOpt(kcnt));aMin2(qOpt(kcnt))]; 
    timeElapsed = toc
end


aOpt

hold on;
plot(kappaSet, qOpt,'o-')
plot(kappaSet, mean(aOpt),'ro-')

save qOptaOptdata.mat qOpt aOpt kappaSet

legend('optimal qO','optimal a')
xlabel('kappa')

totalTime = toc
