clear all;
close all;


load nVsq0simL2Bayes_kappa.5_lambda1.mat

hold on;

plot(nSet, qThyVal*ones(size(nSet)));
error(nSet, mean(qErr),std(qErr));


xlabel('number of data points')
ylabel('Averaged square error in single coefficient')