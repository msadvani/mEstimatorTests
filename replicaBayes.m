clear all;
close all;

lambda = 6;
kappa = .5;



sigma = @(x) lambda*abs(x);


g = @(x)(1/(sqrt(2*pi)))*exp(-x.^2/2); %distribution over beta0


%form a matrix over s0 and q0...
s0 - prox_f(s0 - eta*kappa*(q0 + 1))

constr1 =  @(q0,c) sum()










%%%%%%%%%%%%%%%%%%%%%%%


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

    
    
    
    bMin = fminunc(@(b)(norm(y-X*b,2)^2 + lambda*norm(b,1)),randn(p,1));

    qErr(cnt) = (1/p)*norm(bMin - beta0,2)^2;
      
end


disp('Numerical Experiment: 60% confidence interval')
[mean(qErr)-std(qErr),mean(qErr)+std(qErr)]


disp('Theory')
qThy(kappa)

