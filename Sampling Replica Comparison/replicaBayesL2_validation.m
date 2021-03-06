clear all;
close all;



fNoise = @(x) exp(x)./((1+exp(x)).^2); %noise distribution
%g = @(x) (1/2)*exp(-abs(x)); %coefficient distrib.

g = @(x)probDist(x)

varNoise = sumIntIndef(@(x) x.^2.*fNoise(x), 2,10,.00001)
varCoeff = sumIntIndef(@(x) x.^2.*g(x), 2,10,.00001)



lambda = 1.5;
kappa = .9;


% b = 1-kappa-kappa*lambda;
% c = -kappa*lambda;
%  
% d = (-b+sqrt(b^2-4*c))/2;
% 
% %Note this needs to be changed based on the variance in f and g
% qThy  = @(k) ((d^2 + k*2)./((1+d)^2 - k));




qThyVal = qThyL2(kappa, lambda, varCoeff, varNoise);

%% Simulation
n=800;
p = round(n*kappa);
s = 1;



numSim = 8;
qErr = zeros(1,numSim);
for cnt = 1:numSim
    X = (1/sqrt(p))*randn(n,p); %random unit normal design matrix
    
%     beta0 = randn(p,1);
    sign = 2*round(rand(p,1))-1;
    beta0 = sign.*exprnd(s,p,1);  %double exponential


%Put in laplacian here

%     sign = 2*round(rand(n,1))-1;
%     epsilon = sign.*exprnd(s,n,1);  %double exponential
    
 

    epsilon = random('logistic',0,1,n,1);
        
    %epsilon = randn(n,1);
    
    y = X*beta0 + epsilon;

    
    
    
    bMin = fminunc(@(b)(norm(y-X*b,2)^2 + lambda*norm(b,2)^2),randn(p,1));

    qErr(cnt) = (1/p)*norm(bMin - beta0,2)^2;
      
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


disp('Numerical Experiment: Estimate of 95% confidence interval')
[mean(qErr)-2*std(qErr),mean(qErr)+2*std(qErr)]


disp('Theory')
qThyVal

