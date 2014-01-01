fNoise = @(x) (1/2)*exp(-abs(x));


numKappa = 20;
kappaSet = linspace(.05,.9,numKappa)


ratio = zeros(numKappa,1);
for kcnt = 1:numKappa
    [kcnt, numKappa]
    kappa = kappaSet(kcnt);
    %% Optimal error

    L0 = 10; %initial interval length
    N0=10^3; %initial num point in integral approx
    dblErr = .001; %permissible error in I
    IQ = @(x) IQfunc(fNoise,x,L0,N0, dblErr);

    fToMin = @(z)(IQ(z) - kappa).^2;

    qOpt = gridMinSearch(fToMin,0,5,100,.001);


    qL2 = @(k) 2*k./(1-k);
    
    ratio(kcnt) = qOpt./qL2(kappa);
end


plot(kappaSet, ratio);