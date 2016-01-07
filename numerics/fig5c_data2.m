% synex: scan alpha_k
clear all

% dimensionless parameters
aks = logspace(-2,1,100);
am = 0.3;
bk = 15;
bm = 10;
rho = 0.5;
mu = 1;
h = 2;
p = 2;

% parameter that sets molecule number
kk = 10;

% parameter that sets timescale
lambda = 1;

% dimensionful parameters, for stochastic description
alphaks = aks*kk*lambda;
alpham = am*rho*kk*lambda;
betak = bk*kk*lambda;
betam = bm*rho*kk*lambda;
km = rho*kk;
delta = mu*lambda/rho/kk;

% stochastic parameters
N = 50; nvec = (0:N)'; % ComK molecule number
M = 70; mvec = (0:M)'; % MecA molecule number
I0 = 2000; % number of iterations: first ak value
I1 = 300; % number of iterations: subsequent ak values
Pstart = ones((N+1)*(M+1),1); % starting guess for P

% initialize fixed point structures with null value (-10)
A = length(aks);
Kstar = -10*ones(5,A); % allow at most 5 fixed points
Mstar = -10*ones(5,A);

% loop over ak values
for i = 1:A
  tic
  i
  ak = aks(i);
  alphak = alphaks(i);

  % nullclines
  Kbar = logspace(-3,1,1e4);
  g = ak + bk*Kbar.^h./(1+Kbar.^h);
  q = am + bm*Kbar.^p./(rho^p+Kbar.^p);
  Mbar1 = q;
  Mbar2 = (g./Kbar-1)/mu;

  % fixed points
  DeltaMbar = Mbar1-Mbar2;
  istar = find(diff(sign(DeltaMbar)));
  Kstar(1:length(istar),i) = Kbar(istar)*kk;
  Mstar(1:length(istar),i) = Mbar1(istar)*km;

  % stochastic
  gn = alphak+betak*nvec.^h./(kk^h+nvec.^h);
  qn = alpham+betam*nvec.^p./(km^p+nvec.^p);
  rm = delta*mvec+lambda;
  sbar = lambda;
  %'brute running'
  if i == 1
    [pnm(:,:,i),Pstart,djs0] = ...
        brute_SynEx_start_conv(gn,qn,rm,sbar,N,M,I0,Pstart);
  else
    [pnm(:,:,i),Pstart,djs(:,i)] = ...
        brute_SynEx_start_conv(gn,qn,rm,sbar,N,M,I1,Pstart);
  end
  toc
end

dd = '';
save([dd 'synex_scan.mat'],...
     'aks','N','M','pnm','djs0','djs','Kstar','Mstar')
