% two-Poisson method
clear all

dd = '';
op = optimset('display','final','maxfunevals',1e4,...
              'tolx',1e-8,'tolfun',1e-8,'maxiter',1e4);

% NATIVE
load([dd 'native_scan2.mat']);
aks_native1 = aks;
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
  pi2_native_dkl1(i) = 1-pi1star;

  % plot
  clf
  plot(nvec,pn,'.-',nvec,pnstar,'k')
  xlabel('n (ComK)')
  ylabel('p_n')
  title(['NATIVE: a_k = ' num2str(aks(i))])
  drawnow;% pause
end

load([dd 'native_scan.mat']);
aks_native2 = aks;
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
  pi2_native_dkl2(i) = 1-pi1star;

  % plot
  clf
  plot(nvec,pn,'.-',nvec,pnstar,'k')
  xlabel('n (ComK)')
  ylabel('p_n')
  title(['NATIVE: a_k = ' num2str(aks(i))])
  drawnow;% pause
end

aks_native = [aks_native1 aks_native2];
pi2_native_dkl = [pi2_native_dkl1 pi2_native_dkl2];

% SYNEX
load([dd 'synex_scan2.mat']);
aks_synex1 = aks;
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
  pi2_synex_dkl1(i) = 1-pi1star;

  % plot
  clf
  plot(nvec,pn,'.-',nvec,pnstar,'k')
  xlabel('n (ComK)')
  ylabel('p_n')
  title(['SYNEX: a_k = ' num2str(aks(i))])
  drawnow;% pause
end

load([dd 'synex_scan.mat']);
aks_synex2 = aks;
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
  pi2_synex_dkl2(i) = 1-pi1star;

  % plot
  clf
  plot(nvec,pn,'.-',nvec,pnstar,'k')
  xlabel('n (ComK)')
  ylabel('p_n')
  title(['SYNEX: a_k = ' num2str(aks(i))])
  drawnow;% pause
end

aks_synex = [aks_synex1 aks_synex2];
pi2_synex_dkl = [pi2_synex_dkl1 pi2_synex_dkl2];

% SAVE
save([dd 'pi2_dkl_full.mat'],'pi2_native_dkl','pi2_synex_dkl','aks_native',...
    'aks_synex')