fNoise = @(x) (1/2)*exp(-abs(x));


numKappa = 30;
kappaSet = linspace(.05,.95,numKappa);

 qL2 = @(k) 2*k./(1-k);

qOpt = zeros(numKappa,1);
ratio = zeros(numKappa,1);


for kcnt = 1:numKappa
    [kcnt, numKappa]
    kappa = kappaSet(kcnt);
    %% Optimal error

    L0 = 10; %initial interval length
    N0=10^3; %initial num point in integral approx
    dblErr = .0001; %permissible error in I
    
    I = @(x) calculateInfoMat(fNoise,x,L0,N0, dblErr);
       
    fToMin = @(z)(z.*I(z) - kappa);
    qOpt(kcnt) = findZeroBB(fToMin,.1,5,.001);
       
    ratio(kcnt) = qOpt(kcnt)./qL2(kappa);
end


qOptPlot=plot(kappaSet, qOpt,'o-');
set(qOptPlot, 'LineWidth',1.5)

xlabel('kappa')
ylabel('q0')

save qOptNonBayes.mat qOpt kappaSet

%figure
%plot(kappaSet, ratio,'r-');
