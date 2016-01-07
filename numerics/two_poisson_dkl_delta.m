function [dkl,pfit] = two_poisson_dkl_delta(theta,vec,p,L)

lambda1 = L/(1+exp(-theta(1)));
deltalambda = L/(1+exp(-theta(2)));
lambda2 = lambda1+deltalambda;
pi1 = 1/(1+exp(-theta(3)));

pfit = pi1*poisspdf(vec,lambda1) + (1-pi1)*poisspdf(vec,lambda2);
dkl = sum(p.*log(p./pfit));