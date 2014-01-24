clear all;
close all;

numSim = 10;



qErr = zeros(1,numSim);




%definition of probability distributions
gauss = @(x,sigma) (2*pi*sigma.^2).^(-1/2)*exp(-x.^2./(2*sigma.^2));

sigNoise = 1;
fNoise = @(x) gauss(x,sigNoise);
%fNoise = @(x) 1/2*exp(-abs(x));

g= @(x) 1/2*exp(-abs(x));

varNoise = sumIntIndef(@(x) x.^2.*fNoise(x), 2,10,.00001);
varCoeff = sumIntIndef(@(x) x.^2.*g(x), 2,10,.00001)


%Also try with form 2
form =1;
lambda =1;


kappa = 2.5;
    


q0Replica =compute_q0(kappa, fNoise, g,form, lambda); 



%% Simulation
n=1600;
p = round(n*kappa);
s = 1;

disp('Simulation')


for cnt = 1:numSim
    [cnt,numSim]
    X = (1/sqrt(p))*randn(n,p); %random unit normal design matrix

%     beta0 = randn(p,1);
    sign = 2*round(rand(p,1))-1;
    beta0 = sign.*exprnd(s,p,1);  %double exponential


%Put in laplacian here

%     sign = 2*round(rand(n,1))-1;
%     epsilon = sign.*exprnd(s,n,1);  %double exponential



    %epsilon = random('logistic',0,1,n,1);

    epsilon = randn(n,1)*sigNoise;
    %sign = 2*round(rand(n,1))-1;
    %epsilon = sign.*exprnd(s,n,1);  %double exponential


    y = X*beta0 + epsilon;


    if (form==1)
        bMin = fminunc(@(b)(1/2*norm(y-X*b,2)^2 + lambda*norm(b,1)),randn(p,1));
    elseif (form ==2)
        bMin = fminunc(@(b)((1/2)*norm(y-X*b,2)^2 + (1/2)*lambda*norm(b,2)^2),randn(p,1));
    end


    qErr(cnt) = (1/p)*norm(bMin - beta0,2)^2;
end

disp('Numerical Experiment: Estimate of 65% confidence interval')
[mean(qErr)-std(qErr),mean(qErr)+std(qErr)]

disp('Replica Result')
q0Replica











