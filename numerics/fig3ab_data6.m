% synex: monostable regime
clear all

% dimensionless parameters
ak = 10;
am = 0.3;
bk = 15;
bm = 10;
rho = 0.5;
mu = 1;
h = 2;
p = 2;

% timeseries parameters
X0 = [1e-1; 1e1]; % X = [Kbar; Mbar]
tspan = [0 40];

% timeseries
[tau,Xbar] = ode15s('synex_dim_ode',tspan,X0,[],...
                    ak,am,bk,bm,rho,mu,h,p);

% parameter that sets molecule number (units: molecules)
kk = 10;

% parameter that sets timescale (units: 1/sec)
lambda = 1e-4;

% dimensionful parameters, for stochastic description
alphak = ak*kk*lambda
alpham = am*rho*kk*lambda;
betak = bk*kk*lambda;
betam = bm*rho*kk*lambda;
km = rho*kk;
delta = mu*lambda/rho/kk;

% stochastic parameters
N = 50; nvec = (0:N)'; % ComK molecule number
M = 70; mvec = (0:M)'; % MecA molecule number
I0 = 2000; % number of iterations
Pstart = ones((N+1)*(M+1),1); % starting guess for P

% stationary distribution
gn = alphak+betak*nvec.^h./(kk^h+nvec.^h);
qn = alpham+betam*nvec.^p./(km^p+nvec.^p);
rm = delta*mvec+lambda;
sbar = lambda;
'brute running'
pnm = brute_SynEx_start_conv(gn,qn,rm,sbar,N,M,I0,Pstart);

t = tau/lambda;
K = Xbar(:,1)*kk;
pn = sum(pnm,2);

% save
dd = '';
save([dd 'synex_monostable.mat'],'t','K','pn','nvec')