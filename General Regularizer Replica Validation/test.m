kappa = .3;
fNoise = @(x) GammaDist(x);
g = @(x) .5*exp(-abs(x));


compute_q0(kappa,fNoise, g,1,2)


