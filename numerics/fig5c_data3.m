% two-Poisson method
clear all

dd = '';
op = optimset('display','final','maxfunevals',1e4,...
              'tolx',1e-8,'tolfun',1e-8,'maxiter',1e4);

% NATIVE
load([dd 'native_scan.mat']);
aks_native = aks;
nvec = (0:N)';

for i = 1:length(aks)
  pn = sum(pnm(:,:,i),2);

  % guess
  mu = nvec'*pn;
  lambda1 = mu-rand;
  deltalambda = 2*rand;
  pi1 = rand;
  theta0 = [log(lambda1/(N-lambda1)) ...
            log(deltalambda/(N-deltalambda)) ...
            log(pi1/(1-pi1))];
  
  % fit
  thetastar = fminsearch('two_poisson_dkl_delta',theta0,op,...
                         nvec,pn,N);
  [dklstar,pnstar] = two_poisson_dkl_delta(thetastar,nvec,pn,N);
  lambda1star = N/(1+exp(-thetastar(1)));
  deltalambdastar = N/(1+exp(-thetastar(2)));
  lambda2star = lambda1star+deltalambdastar;
  pi1star = 1/(1+exp(-thetastar(3)));
  pi2_native_dkl(i) = 1-pi1star;

  % plot
  clf
  plot(nvec,pn,'.-',nvec,pnstar,'k')
  xlabel('n (ComK)')
  ylabel('p_n')
  title(['NATIVE: a_k = ' num2str(aks(i))])
  drawnow;% pause
end

% SYNEX
load([dd 'synex_scan.mat']);
aks_synex = aks;
nvec = (0:N)';

for i = 1:length(aks)
  pn = sum(pnm(:,:,i),2);

  % guess
  mu = nvec'*pn;
  lambda1 = mu-rand;
  deltalambda = 2*rand;
  pi1 = rand;
  theta0 = [log(lambda1/(N-lambda1)) ...
            log(deltalambda/(N-deltalambda)) ...
            log(pi1/(1-pi1))];

  % fit
  thetastar = fminsearch('two_poisson_dkl_delta',theta0,op,...
                         nvec,pn,N);
  [dklstar,pnstar] = two_poisson_dkl_delta(thetastar,nvec,pn,N);
  lambda1star = N/(1+exp(-thetastar(1)));
  deltalambdastar = N/(1+exp(-thetastar(2)));
  lambda2star = lambda1star+deltalambdastar;
  pi1star = 1/(1+exp(-thetastar(3)));
  pi2_synex_dkl(i) = 1-pi1star;

  % plot
  clf
  plot(nvec,pn,'.-',nvec,pnstar,'k')
  xlabel('n (ComK)')
  ylabel('p_n')
  title(['SYNEX: a_k = ' num2str(aks(i))])
  drawnow;% pause
end

% SAVE
save([dd 'pi2_dkl.mat'],'pi2_native_dkl','pi2_synex_dkl')