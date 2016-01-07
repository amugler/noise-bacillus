% native: low induction
% computes f for experimental distributions using 2-poisson method
clear all

plotit = 1;
dd = '';
fd = '';
load([dd 'pixel_data.mat']);
load([dd 'native_fit_cdf.mat']);

op = optimset('display','final','maxfunevals',1e4,...
              'tolx',1e-8,'tolfun',1e-8,'maxiter',1e4);

% dkl with 2 poissons to get weights
for i = 1:4
  x = double(native.px{i})';

  % experimental PDF with calibration
  if i == 1
    x0 = mode(x);
  end
  x_ = x-x0;
  x_(find(x_ < 0)) = 0;
  n_ = x_/Xstar;
  N_ = ceil(max(n_));
  n_vec = (0:N_)';
  pn_ = hist(n_,n_vec);
  pn_ = pn_'/sum(pn_);

  % guess
  mu = n_vec'*pn_;
  lambda1 = mu-.3*rand;
  deltalambda = 2*rand;
  pi1 = rand;
  theta0 = [log(lambda1/(N_-lambda1)) ...
            log(deltalambda/(N_-deltalambda)) ...
            log(pi1/(1-pi1))];

  % fit
  thetastar = fminsearch('two_poisson_dkl_delta_plot',theta0,op,...
                         n_vec,pn_,N_,plotit);
  [dklstar,pn_star] = ...
      two_poisson_dkl_delta_plot(thetastar,n_vec,pn_,N_,plotit);
  lambda1star = N_/(1+exp(-thetastar(1)));
  deltalambdastar = N_/(1+exp(-thetastar(2)));
  lambda2star = lambda1star+deltalambdastar;
  pi1star = 1/(1+exp(-thetastar(3)));
  pi2_native_dkl_exp(i) = 1-pi1star;
end

% save
save([dd 'pi2_native_dkl_exp_cdf.mat'],'pi2_native_dkl_exp')
