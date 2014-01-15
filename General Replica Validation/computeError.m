function [qOut, cOut] = computeError(kappa, fNoise, g,tol )
%Computes q0 and c given kappa, noise, coefficient distrib, and tolerance

%These paramters may have to be tweaked...they are the grid spacing on the
%integrals being approximated.

N=5000;
L=13;

%gauss = @(x,sigma) (2*pi*sigma.^2).^(-1/2)*exp(-x.^2./(2*sigma.^2));

varNoise = sumIntIndef(@(x) x.^2.*fNoise(x), 2,10,.00001);
%varCoeff = sumIntIndef(@(x) x.^2.*g(x), 2,10,.00001)

c0Hat = @(q0,c) (1/(2*kappa*(1+c).^2))*(q0 + varNoise);
cHat = @(q0,c) (1./(2*kappa*(1+c)));
a = @(q0) (q0 + varNoise)*kappa;




%%%[q0Thy,cThy] = qThyL2(kappa, 1, varCoeff, varNoise);

max_pErr = .01; 
F1 = @(q0,c)computeConstraint1(q0,c,fNoise, g,kappa,N,L,max_pErr);
F2 = @(q0,c)computeConstraint2(q0,c,fNoise, g,kappa,N,L,max_pErr);



%lambda =1;
%F1Thy = @(q0,c) -q0 + (lambda^2*kappa^2*(1+c).^2*varCoeff + a(q0))./(1+lambda*kappa*(1+c)).^2
%F2Thy = @(q0,c) c./(1+c) - kappa./(1+lambda*kappa*(1+c))


c1Min = @(q) findZeroBB_cutoff(@(c)F1(q,c),.01,2,tol);
c2Min = @(q) findZeroBB_cutoff(@(c)F2(q,c),.01,2,tol);


qOut = findZeroBB(@(q) (c1Min(q)-c2Min(q)),.01,2,tol);

cOut = [c1Min(qOut),c2Min(qOut)];

end

