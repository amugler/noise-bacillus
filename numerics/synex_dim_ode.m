function F = synex_dim_ode(t,X,ig,ak,am,bk,bm,rho,mu,n,p)

Kbar = X(1);
Mbar = X(2);
g = ak + bk*Kbar^n/(1+Kbar^n);
q = am + bm*Kbar^p/(rho^p+Kbar^p);
r = mu*Mbar+1;
s = 1;
F(1,1) = g - r*Kbar;
F(2,1) = q - s*Mbar;
