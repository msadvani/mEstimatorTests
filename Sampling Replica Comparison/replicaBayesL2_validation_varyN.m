clear all;
close all;



lambda = 1;
kappa = .5;


b = 1-kappa-kappa*lambda;
c = -kappa*lambda;
d = (-b+sqrt(b^2-4*c))/2;

 
qThyVal = (d^2 + kappa*2)./((1+d)^2 - kappa);

numSim = 10; %number of simulations for each value of n


nSet = [100:100:2000];
qErr = zeros(numSim, length(nSet));

for nCnt = 1:length(nSet)
    [nCnt,length(nSet)];
    n=nSet(nCnt);
    p = round(n*kappa);
    s = 1;



  
    
    for simCnt = 1:numSim
        X = (1/sqrt(p))*randn(n,p); %random unit normal design matrix
        beta0 = randn(p,1);

        sign = 2*round(rand(n,1))-1;
        epsilon = sign.*exprnd(s,n,1);  %double exponential
        %epsilon = randn(n,1);

        y = X*beta0 + epsilon;

        bMin = fminunc(@(b)(norm(y-X*b,2)^2 + lambda*norm(b,2)^2),randn(p,1));

        qErr(simCnt,nCnt) = (1/p)*norm(bMin - beta0,2)^2;

    %       rho = @(x)(1/2)*x.^2;
    %       
    %       %sigma = @(x)zeros(size(x)); %Easy case
    %       sigma = @(x) (1/2)*lambda*(x.^2);
    %       
    %       tic
    %       [bMinS]=fminunc(@(b)sum(rho(y-X*b)), randn(p,1));
    %       %[bMinS] = fminunc(@(b) sum(rho(y-X*b)) + sum(sigma(b)),randn(p,1));
    %       toc
    %       
    %       rhoErr(cnt) = (1/sqrt(P))*norm(bMinS - beta0,2);
    end
end

%disp('Numerical Experiment: Estimate of 95% confidence interval')
%[mean(qErr)-2*std(qErr),mean(qErr)+2*std(qErr)]

%disp('Theory')
%qThy(kappa)

hold on;
plot(nSet, qThyVal*ones(size(nSet)))
plot(nSet, mean(qErr)-std(qErr),'r--')
plot(nSet, mean(qErr)+std(qErr),'r--')

legend('prediction','1 std in simulation')

