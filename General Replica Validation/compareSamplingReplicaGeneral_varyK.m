clear all;
close all;

numSim = 10;

numK=16;
kappaSet = linspace(.1,1.5,numK);
q0Replica = zeros(size(kappaSet));

qErr = zeros(numSim,numK);


tic
for kcnt=1:numK
    disp('Big Loop')
    [kcnt,numK]
    kappa = kappaSet(kcnt);
    
    gauss = @(x,sigma) (2*pi*sigma.^2).^(-1/2)*exp(-x.^2./(2*sigma.^2));

    sigNoise = 1;
    fNoise = @(x) gauss(x,sigNoise);
    g = @(x) 1/2*exp(-abs(x));

    varNoise = sumIntIndef(@(x) x.^2.*fNoise(x), 2,10,.00001);
    varCoeff = sumIntIndef(@(x) x.^2.*g(x), 2,10,.00001)

    c0Hat = @(q0,c) (1/(2*kappa*(1+c).^2))*(q0 + varNoise);
    cHat = @(q0,c) (1./(2*kappa*(1+c)));

    a = @(q0) (q0 + varNoise)*kappa;


    N=5000;
    L=13;


    %[q0Thy,cThy] = qThyL2(kappa, 1, varCoeff, varNoise)
    %Now we need to compare this with the analytic solution already worked out!

    %computeConstraint1(1,1,fNoise, g,kappa,N,L,.1)
    F1 = @(q0,c)computeConstraint1(q0,c,fNoise, g,kappa,N,L,.01)
    F2 = @(q0,c)computeConstraint2(q0,c,fNoise, g,kappa,N,L,.01)


    % lambda =1;
    % F1Thy = @(q0,c) -q0 + (lambda^2*kappa^2*(1+c).^2*varCoeff + a(q0))./(1+lambda*kappa*(1+c)).^2
    % 
    % F2Thy = @(q0,c) c./(1+c) - kappa./(1+lambda*kappa*(1+c))


    Ftot = @(q0,c) (F1(q0,c)^2+F2(q0,c)^2);

    [q0,c,err]=gridMin2d(Ftot,0,0,3,3,.01);
    
    q0Replica(kcnt) =q0; 




    %% Simulation
    n=800;
    p = round(n*kappa);
    s = 1;
    lambda =1;

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

        y = X*beta0 + epsilon;




        %bMin = fminunc(@(b)(norm(y-X*b,2)^2 + lambda*norm(b,2)^2),randn(p,1));

        bMin = fminunc(@(b)(1/2*norm(y-X*b,2)^2 + lambda*norm(b,1)),randn(p,1));


        qErr(cnt,kcnt) = (1/p)*norm(bMin - beta0,2)^2;
    end

%     disp('Numerical Experiment: Estimate of 65% confidence interval')
%     [mean(qErr)-std(qErr),mean(qErr)+std(qErr)]
% 
%     disp('Replica Result')
%     q0

end


save generalValidateReplica.mat kappaSet qReplica qErr

hold on;
kappaMat = repmat(kappaSet, numSim, 1);
kappaSetVec = reshape(kappaMat,1,numKappa*numSim);
qErrVec = reshape(qErr,1,numKappa*numSim);


scatter(kappaSetVec, qErrVec,'r.');

pThy = plot(kappaSet, qReplica,'k');

set(pThy,'Linewidth',2.5);

legend('Theory','Simulations')

xlabel('kappa');
ylabel('error')
toc




title('Validating results with Gaussian Noise and L1 regularizer')




















