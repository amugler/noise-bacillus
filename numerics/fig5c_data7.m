% native: high induction
% fits experimental pixel CDFs with stochastic model
% optimizes over one calibration parameter:
%   X = pixel intensity of one molecule
% cost function is sum-of-squares per data point,
%   minimized over a_k values for each IPTG level,
%   then summed over all IPTG levels
clear all

op = optimset('display','iter','maxfunevals',1e4,...
              'tolx',1e-4,'tolfun',1e-4,'maxiter',1e4);
dd = '';

% boundaries for calibration parameter
Xmin = 1e0;
Xmax = 2e2;

% optimize
Xstar = fminbnd('native_cost3',Xmin,Xmax,op);
[Ssum_min,kstar,S,n1,n2,c1,c2] = native_cost3(Xstar);

% save
save([dd 'native_fit_cdf3.mat'],'Xstar','kstar','S','n1','n2','c1','c2')
