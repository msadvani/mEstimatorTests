fNoise = @(x) (1/2)*exp(-abs(x));


numKappa = 50;
kappaSet = linspace(.02,.98,numKappa);


ratio = zeros(numKappa,1);
for kcnt = 1:numKappa
    [kcnt, numKappa]
    kappa = kappaSet(kcnt);
    %% Optimal error

    L0 = 10; %initial interval length
    N0=10^3; %initial num point in integral approx
    dblErr = .0001; %permissible error in I
    IQ = @(x) IQfunc(fNoise,x,L0,N0, dblErr);
    
    fToMin = @(z)(IQ(z) - kappa);
    qOpt = findZeroBB(fToMin,.1,5,.0001);

%     fToMin = @(z)(IQ(z) - kappa).^2;
% 
%     qOpt = gridMinSearch(fToMin,0,5,100,.001);


    qL2 = @(k) 2*k./(1-k);
    
    ratio(kcnt) = qOpt./qL2(kappa);
end


plot(kappaSet, ratio,'*');