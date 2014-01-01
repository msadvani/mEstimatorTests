clear all;
close all;

kappa = .5;
lambda = 1;


fNoise = @(x) (2*pi)^(-1/2)*exp(-x.^2/2);
g = @(x) (1/2)*exp(-abs(x));


%% Theory for Optimal error

L0 = 10; %initial interval length
N0=10^3; %initial num point in integral approx
dblErr = .0001; %permissible error in I


%Ooops...maybe need to make the vector veions more general...
I = @(x) calculateInfoMat(fNoise,x,L0,N0, dblErr);
J = @(x) calculateInfoMat(g,x, L0, N0, dblErr);



%plotting the function
%a = linspace(0,50,200);
%plot(a,a - a.^2.*J(a))

%two constraints
F1 = @(q,a) (a - a.^2.*J(a) - q).^2;
%F2 = @(q,a) (I(q).*a - kappa).^2;


numA = 50;
numQ = 50;
aSet =repmat(linspace(0,400,numA),numQ,1);
qSet = repmat(linspace(0,10,numQ)',1,numA);

imagesc(F1(qSet, aSet))
colorbar

%plot(aSet,F1(qSet(4),aSet))

% for cnt =1:numQ
%    F1(qSet(cnt),aSet))
% end





%numR =100;
%numC = 200;

%aSet = repmat(linspace(0,5,numR)',1,numC);
%qSet = repmat(linspace(0,5,numC),numR,1);


%Fmat = F1(qSet,aSet)





%fToMin = @(z)(IQ(z) - kappa).^2;

%qOpt = gridMinSearch(fToMin,0,5,100,.005);








% %% Simulation
% n=500;
% p = round(n*kappa);
% s = 1;
% 
% 
% 
% 
% numSim = 5;
% qErr = zeros(1,numSim);
% for cnt = 1:numSim
%     X = (1/sqrt(p))*randn(n,p); %random unit normal design matrix
%     beta0 = randn(p,1);
% 
%     sign = 2*round(rand(n,1))-1;
%     epsilon = sign.*exprnd(s,n,1);  %double exponential
%     %epsilon = randn(n,1);
%     
%     y = X*beta0 + epsilon;
%     
%     bMin = fminunc(@(b)(norm(y-X*b,2)^2 + lambda*norm(b,1)),randn(p,1));
% 
%     qErr(cnt) = (1/p)*norm(bMin - beta0,2)^2;
% end
% 
% disp('Numerical Experiment: 60% confidence interval')
% [mean(qErr)-std(qErr),mean(qErr)+std(qErr)]
% 
% 
% disp('Theory')
% qThy(kappa)
% 
