function [q0_out] = compute_q0(kappa, fNoise, g,form,lambda)
%Computes an estimate of analytic error of q0, form = type of convex
%function, lambda = regularization parameter

varNoise = sumIntIndef(@(x) x.^2.*fNoise(x), 5,100,.00001);

N=5000;
L=13;

%form = 1; %prox of what form? (1)- L1 norm
%lambda =1; %regularization paramter


%computeConstraint1(1,1,fNoise, g,kappa,N,L,.1)
F1 = @(q0,c)computeConstraint1(q0,c,fNoise, g,kappa,form, lambda, N,L,.01);
F2 = @(q0,c)computeConstraint2(q0,c,fNoise, g,kappa,form, lambda, N,L,.01);

Ftot = @(q0,c) (F1(q0,c)^2+F2(q0,c)^2);

[q0,c,err]=gridMin2d(Ftot,0,0,3,3,.01)


q0_out = q0;

end

