% native: scan alpha_k (lower values)
clear all

% dimensionless parameters
aks = logspace(-5,-4,30);
as = 0.0;
bk = 0.3;
bs = 3;
k0 = 0.2;
k1 = 1/30;
Deltak = 0.1;
Deltas = 0.1;
h = 2;
p = 5;

% parameters that set molecule number
Gammak = 100;
Gammas = 1;

% parameters that set timescale
deltak = 1;
deltas = 1;

% dimensionful parameters, for stochastic description
alphaks = aks*Gammak*deltak;
alphas = as*Gammas*deltas;
betak = bk*Gammak*deltak;
betas = bs*Gammas*deltas;
lambdak = Deltak*deltak;
lambdas = Deltas*deltas;
kk = k0*Gammak;
ks = k1*Gammak;

% stochastic parameters
N = 80; nvec = (0:N)'; % ComK molecule number
M = 50; mvec = (0:M)'; % ComS molecule number
nmat = nvec*ones(1,M+1);
mmat = ones(N+1,1)*mvec';
I0 = 2000; % number of iterations: first ak value
I1 = 300; % number of iterations: subsequent ak values
Pstart = ones((N+1)*(M+1),1); % starting guess for P

% initialize fixed point structures with null value (-10)
A = length(aks);
Kstar = -10*ones(8,A); % allow at most 8 fixed points
Sstar = -10*ones(8,A);

% loop over ak values
for i = 1:A
  i
  ak = aks(i);
  alphak = alphaks(i);

  % nullclines
  Kbar = logspace(-3,1,1e6);
  g = ak + bk*Kbar.^h./(k0^h+Kbar.^h);
  q = as + bs./(1+(Kbar/k1).^p);
  Sbar1 = 1./(g./Kbar-Deltak)-Kbar-1;
  B = 1+(1+Kbar)*Deltas-q;
  C = (1+Kbar).*q;
  Sbar2 = (sqrt(B.^2+4*Deltas*C)-B)/2/Deltas;

  % fixed points
  DeltaSbar = Sbar1-Sbar2;
  istar = find(diff(sign(DeltaSbar)));
  Kstar(1:length(istar),i) = Kbar(istar)*Gammak;
  Sstar(1:length(istar),i) = Sbar1(istar)*Gammas;

  % stochastic
  gn = alphak+betak*nvec.^h./(kk^h+nvec.^h);
  qn = alphas+betas./(1+(nvec/ks).^p);
  rnm = deltak./(1+nmat/Gammak+mmat/Gammas)+lambdak;
  snm = deltas./(1+nmat/Gammak+mmat/Gammas)+lambdas;
  %'brute running'
  if i == 1
    [pnm(:,:,i),Pstart,djs0] = ...
        brute_Native_start_conv(gn,qn,rnm,snm,N,M,I0,Pstart);
  else
    [pnm(:,:,i),Pstart,djs(:,i)] = ...
        brute_Native_start_conv(gn,qn,rnm,snm,N,M,I1,Pstart);
  end

end

dd = '';
save([dd 'native_scan2.mat'],...
     'aks','N','M','pnm','djs0','djs','Kstar','Sstar')
