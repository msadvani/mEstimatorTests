close all;
clear all;

load OptM_LapPred.mat
load replicaMapError_lambda=1.mat
load replicaL2Error_lambda=1.mat



hold on;


p1=plot(kappaSet, q0ReplicaMAP,'b');
p2= plot(kappaSet, q0ReplicaL2,'r');
p3 = plot(kappaSet, qOpt,'k--');

set(p1,'Linewidth',2)
set(p2,'Linewidth',2)
set(p3,'Linewidth',2);

legend('MAP','L2 Regularizer','Optimal M-estimator');

%xlabel('\kappa  (P/N)')
%ylabel('q0 MSE per predictor')

%title('Comparison of Asymptotic performance of MAP, L2, and Optimal M-Estimator Error')

%Note you can add in the data measurements if you like, but it may make the
%plot harder to read.


tlhand = get(gca,'title')
set(tlhand,'string','Asymptotic MSE of MAP, L2 regularizer, and Optimal M-Estimator','fontsize',18)

xlhand = get(gca,'xlabel')
set(xlhand,'string','\kappa (P/N)','fontsize',16) 


ylhand = get(gca,'ylabel')
set(ylhand,'string','Mean Squared Error per Coefficient (MSE)','fontsize',16) 

%Set font size of numbering on axes
set(gca,'FontSize',14)












%% Fix and use for scatter plots
% numSim = size(qErrMAP,1);
% numK = size(qErrMAP,2)
% kappaMat = repmat(kappaSet, numSim, 1);
% kappaSetVec = reshape(kappaMat,1,numK*numSim);
% qErrVecMAP = reshape(qErrMAP,1,numK*numSim);
% 
% scatter(kappaSetVec, qErrVecMAP,'r.');
 
%pThy = plot(kappaSet, q0Replica,'k');
% 
% set(pThy,'Linewidth',2.5);
% 
% %legend('Simulations','Theory')
% 
% xlabel('kappa');
% ylabel('error')
% toc
% 
% 
% 
% 
% title('Validating results with Gaussian Noise Laplacian coefficients and and L1 regularizer, lambda =2')
