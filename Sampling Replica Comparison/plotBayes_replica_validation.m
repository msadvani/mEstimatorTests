clear all;
close all;


load nVsq0simL2Bayes_lapNoise_gaussCoeff_kappa.5_lambda1.mat

hold on;

plot(nSet, qThyVal*ones(size(nSet)));
errorbar(nSet, mean(qErr),std(qErr),'-.or');


xlabel('number of data points')
ylabel('Averaged square error in single coefficient')

title('Comparison of Replica Theory with Numerical Simulations for Estimator p(x) = x^2/2  \sigma(x)=|x| with kappa =.5')

legend('Replica Theory','Numerical Simulations (\pm 1 std)')
