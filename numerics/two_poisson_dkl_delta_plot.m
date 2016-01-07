function [dkl,pfit] = two_poisson_dkl_delta(theta,vec,p,L,plotit)

lambda1 = L/(1+exp(-theta(1)));
deltalambda = L/(1+exp(-theta(2)));
lambda2 = lambda1+deltalambda;
pi1 = 1/(1+exp(-theta(3)));

pfit = pi1*poisspdf(vec,lambda1) + (1-pi1)*poisspdf(vec,lambda2);
nz = find(p);
dkl = sum(p(nz).*log(p(nz)./pfit(nz)));

if plotit
  figure(1); clf; hold on
  bar(vec,p);
  plot(vec,pfit,'m.-')
  drawnow
end
