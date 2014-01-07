%This code compares the optimal L2 based loss function error with the
%predicted optimal erro based on my research. qOpt should always be lower,
%apart from numerical inaccuracy.

clear all;
close all;

kappaSet = [.3,1.3,2.3];
numKappa = length(kappaSet);

%numKappa = 20;
%kappaSet = linspace(.1,2,numKappa);


fNoise = @(x) (2*pi)^(-1/2)*exp(-x.^2/2); %noise 
g = @(x) (1/2)*exp(-abs(x)); %coefficient distrib.

%g = @(x) (2*pi)^(-1/2)*exp(-x.^2/2); %coefficient distrib.


qOpt = zeros(1, numKappa);

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
    F1 = @(q,a) (a - a.^2.*J(a) - q).^2;
    F2 = @(q,a) (I(q).*a - kappa).^2;


    aMin1 = @(q)gridMinSearch(@(x)F1(q,x),0,5,10,.01);
    aMin2 = @(q)gridMinSearch(@(x)F2(q,x),0,5,10,.01);


    d = @(q)(aMin1(q)-aMin2(q)).^2;

    qOpt(kcnt) = gridMinSearchNonVec(d,.05,2,10,.01);

end


%% Compute the integral estimating the error for lasso


%% Theory for a suboptimal case (L2 norm) [This comes from replica equations]

qL2Min = zeros(1,numKappa);
for kcnt = 1:numKappa
    kappa = kappaSet(kcnt);
    f=@(lam)qThyL2(kappa,lam);
    lamMin = gridMinSearch(f,0,10,100,.01);
    qL2Min(kcnt) = f(lamMin);
end

save qOptVsL2data.mat kappaSet qL2Min qOpt


hold on;
plot(kappaSet, qL2Min,'r')
plot(kappaSet, qOpt)


