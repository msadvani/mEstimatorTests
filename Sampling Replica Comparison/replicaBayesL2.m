clear all;
close all;



lambda = 1.5;
kappa = .9;


b = 1-kappa-kappa*lambda;
c = -kappa*lambda;
 
d = (-b+sqrt(b^2-4*c))/2;

 

%Note this needs to be changed based on the variance in f and g
qThy  = @(k) ((d^2 + k*2)./((1+d)^2 - k));


%% Simulation
n=500;
p = round(n*kappa);
s = 1;




numSim = 5;
qErr = zeros(1,numSim);
for cnt = 1:numSim
    X = (1/sqrt(p))*randn(n,p); %random unit normal design matrix
    beta0 = randn(p,1);

    sign = 2*round(rand(n,1))-1;
    epsilon = sign.*exprnd(s,n,1);  %double exponential
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
qThy(kappa)

