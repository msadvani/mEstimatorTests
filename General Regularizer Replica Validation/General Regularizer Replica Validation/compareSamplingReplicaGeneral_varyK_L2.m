clear all;
close all;

numSim = 10;

numK=30;
kappaSet = linspace(.05,4.0,numK);
q0Replica = zeros(size(kappaSet));

qErr = zeros(numSim,numK);




%definition of probability distributions
gauss = @(x,sigma) (2*pi*sigma.^2).^(-1/2)*exp(-x.^2./(2*sigma.^2));

sigNoise = 1;
fNoise = @(x) gauss(x,sigNoise);
%fNoise = @(x) 1/2*exp(-abs(x));

g= @(x) 1/2*exp(-abs(x));

varNoise = sumIntIndef(@(x) x.^2.*fNoise(x), 2,10,.00001);
varCoeff = sumIntIndef(@(x) x.^2.*g(x), 2,10,.00001)


%Also try with form 2
form =2;
lambda =1;


tic
for kcnt=1:numK
    disp('Big Loop')
    [kcnt,numK]
    kappa = kappaSet(kcnt);
    
    
    
    q0Replica(kcnt) =compute_q0(kappa, fNoise, g,form, lambda); 



    %% Simulation
    n=800;
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

        
        qErr(cnt,kcnt) = (1/p)*norm(bMin - beta0,2)^2;
    end

%     disp('Numerical Experiment: Estimate of 65% confidence interval')
%     [mean(qErr)-std(qErr),mean(qErr)+std(qErr)]
% 
%     disp('Replica Result')
%     q0

end

qErrL2 = qErr;
q0ReplicaL2 = q0Replica;
save replicaL2Error_lambda=1.mat kappaSet q0ReplicaL2 qErrL2

hold on;
kappaMat = repmat(kappaSet, numSim, 1);
kappaSetVec = reshape(kappaMat,1,numK*numSim);
qErrVec = reshape(qErr,1,numK*numSim);


scatter(kappaSetVec, qErrVec,'r.');

pThy = plot(kappaSet, q0Replica,'k');

set(pThy,'Linewidth',2.5);

%legend('Simulations','Theory')

xlabel('kappa');
ylabel('error')
toc




title('Validating results with Gaussian Noise Laplacian coefficients and and L2 regularizer, lambda =1')




















