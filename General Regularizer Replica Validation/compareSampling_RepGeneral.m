clear all;
close all;

kappa = .3;
lambda =3;

gauss = @(x,sigma) (2*pi*sigma.^2).^(-1/2)*exp(-x.^2./(2*sigma.^2));

%sigNoise = 1.5;
%fNoise = @(x) gauss(x,sigNoise);
fNoise = @(x) GammaDist(x);

g = @(x) 1/2*exp(-abs(x));


%%%Just a test of function
q0 = compute_q0(kappa,fNoise, g,1,lambda)



varNoise = sumIntIndef(@(x) x.^2.*fNoise(x), 5,100,.00001)

varCoeff = sumIntIndef(@(x) x.^2.*g(x), 2,10,.00001)

c0Hat = @(q0,c) (1/(2*kappa*(1+c).^2))*(q0 + varNoise);
cHat = @(q0,c) (1./(2*kappa*(1+c)));

a = @(q0) (q0 + varNoise)*kappa;


N=5000;
L=13;






%% Simulation
n=800;
p = round(n*kappa);
s = 1;


disp('Simulations')

numSim = 15;
qErr = zeros(1,numSim);
for cnt = 1:numSim
    [cnt,numSim]
    X = (1/sqrt(p))*randn(n,p); %random unit normal design matrix
    
%     beta0 = randn(p,1);
    sign = 2*round(rand(p,1))-1;
    beta0 = sign.*exprnd(s,p,1);  %double exponential


%Put in laplacian here

%     sign = 2*round(rand(n,1))-1;
%     epsilon = sign.*exprnd(s,n,1);  %double exponential
    epsilon = gamrnd(2,1,[n,1]);
 

    %epsilon = random('logistic',0,1,n,1);   
    %epsilon = randn(n,1)*sigNoise;
    
    
    
    y = X*beta0 + epsilon;

    
    %bMin = fminunc(@(b)(norm(y-X*b,2)^2 + lambda*norm(b,2)^2),randn(p,1));

    bMin = fminunc(@(b)(1/2*norm(y-X*b,2)^2 + lambda*norm(b,1)),randn(p,1));
    
    
    qErr(cnt) = (1/p)*norm(bMin - beta0,2)^2;
      
end


disp('Numerical Experiment: Estimate of 65% confidence interval')
[mean(qErr)-std(qErr),mean(qErr)+std(qErr)]

disp('Replica Result')
q0

























