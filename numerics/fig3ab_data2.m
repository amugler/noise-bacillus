% native: oscillatory regime
clear all

% dimensionless parameters
ak = 3.2e-3;
as = 0.0;
bk = 0.3;
bs = 3;
k0 = 0.2;
k1 = 1/30;
Deltak = 0.1;
Deltas = 0.1;
h = 2;
p = 5;

% timeseries parameters
X0 = [3e-2; 5e0]; % X = [Kbar; Sbar]
tspan = [0 200];

% timeseries
[tau,Xbar] = ode15s('native_dim_ode',tspan,X0,[],...
    ak,as,bk,bs,k0,k1,Deltak,Deltas,h,p);

% parameters that set molecule number
Gammak = 100;
Gammas = 1;

% parameters that set timescale (units: 1/sec)
deltak = .001;
deltas = .001;

% dimensionful parameters, for stochastic description
alphak = ak*Gammak*deltak
alphas = as*Gammas*deltas;
betak = bk*Gammak*deltak;
betas = bs*Gammas*deltas;
lambdak = Deltak*deltak;
lambdas = Deltas*deltas;
kk = k0*Gammak;
ks = k1*Gammak;

% stochastic parameters
N = 120; nvec = (0:N)'; % ComK molecule number
M = 50; mvec = (0:M)'; % ComS molecule number
nmat = nvec*ones(1,M+1);
mmat = ones(N+1,1)*mvec';
I0 = 2000; % number of iterations: first ak value
Pstart = ones((N+1)*(M+1),1); % starting guess for P

% stationary distribution
gn = alphak+betak*nvec.^h./(kk^h+nvec.^h);
qn = alphas+betas./(1+(nvec/ks).^p);
rnm = deltak./(1+nmat/Gammak+mmat/Gammas)+lambdak;
snm = deltas./(1+nmat/Gammak+mmat/Gammas)+lambdas;
'brute running'
pnm = brute_Native_start_conv(gn,qn,rnm,snm,N,M,I0,Pstart);

t = tau/deltak;
K = Xbar(:,1)*Gammak;
pn = sum(pnm,2);

% find mode separation by inflection method
% compute second derivative
dpn = diff([0;pn]);
ddpn = diff([0;dpn]);
% smooth out alternations (NOT UNDERSTOOD)
ddpn_ = (ddpn + [0;ddpn(1:end-1)])/2;
% find first four inflection points (btw max, "min", and "max")
ind = find(diff(sign(ddpn_)) ~= 0);
% separator is average of second and thir inflection point
i0 = ceil((ind(2)+ind(3))/2);
f = 1-sum(pn(1:i0))

% save
dd = '';
save([dd 'native_oscillatory.mat'],'t','K','pn','nvec','i0')