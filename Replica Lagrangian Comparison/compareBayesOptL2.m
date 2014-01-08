%This code compares the optimal L2 based loss function error with the
%predicted optimal erro based on my research. qOpt should always be lower,
%apart from numerical inaccuracy.

clear all;
close all;

%kappaSet = [.3,1.3,2.3];
%numKappa = length(kappaSet);
 
numKappa = 25;
kappaSet = linspace(.1,5,numKappa);


%fNoise = @(x) (2*pi)^(-1/2)*exp(-x.^2/2); %noise
%fNoise = @(x) (1/2)*exp(-abs(x)); %noise
fNoise = @(x) probDist(x);



%g = @(x) (2*pi)^(-1/2)*exp(-x.^2/2); %coefficient distrib.
g = @(x) (1/2)*exp(-abs(x)); %coefficient distrib.

varNoise = sumIntIndef(@(x) x.^2.*fNoise(x), 2,10,.00001)
varCoeff = sumIntIndef(@(x) x.^2.*g(x), 2,10,.00001)




qOpt = zeros(1, numKappa);
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
    
    timeElapsed = toc
end


%% Compute the integral estimating the error for lasso


%% Theory for a suboptimal case (L2 norm) [This comes from replica equations]

qL2Min = zeros(1,numKappa);
for kcnt = 1:numKappa
    kappa = kappaSet(kcnt);
    f=@(lam)qThyL2(kappa,lam,varCoeff, varNoise);
    lamMin = gridMinSearch(f,0,10,100,.005);
    qL2Min(kcnt) = f(lamMin);
end

save qOptVsL2data.mat kappaSet qL2Min qOpt


hold on;
plot(kappaSet, qL2Min,'r')
plot(kappaSet, qOpt)

legend('qL2Min','qOpt')
xlabel('kappa')
ylabel('q')

totalTime = toc
