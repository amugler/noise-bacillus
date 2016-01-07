function F = native_dim_ode(t,X,ig,ak,as,bk,bs,k0,k1,Deltak,Deltas,n,p)

Kbar = X(1);
Sbar = X(2);
g = ak + bk*Kbar^n/(k0^n+Kbar^n);
q = as + bs/(1+(Kbar/k1)^p);
r = 1/(1+Kbar+Sbar)+Deltak;
s = 1/(1+Kbar+Sbar)+Deltas;
F(1,1) = g - r*Kbar;
F(2,1) = q - s*Sbar;
