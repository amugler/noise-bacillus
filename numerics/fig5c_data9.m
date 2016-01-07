% synex: intermediate induction
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
Xmin = 1e-1;
Xmax = 1e2;

% optimize
Xstar = fminbnd('synex_cost2',Xmin,Xmax,op);
[Ssum_min,kstar,S,n1,n2,c1,c2] = synex_cost2(Xstar);

% save
save([dd 'synex_fit_cdf2.mat'],'Xstar','kstar','S','n1','n2','c1','c2')
