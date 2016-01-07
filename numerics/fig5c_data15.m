% native: low induction
% computes error bar on a_k for experimental data
clear all

dd = '';
load([dd 'native_scan.mat']);
load([dd 'native_fit_cdf.mat']);

f = .25;
for i = 1:4
  [Smin,kstar] = min(S(i,:));
  Scross(i,1) = Smin*(1+f);
  dS = S(i,:)-Scross(i);
  kcross{i} = find(diff(sign(dS)));
end

figure(1); clf;
semilogx(aks,S,'.-',...
         [min(aks) max(aks)],[Scross Scross],'-')

save([dd 'error_native.mat'],'kcross')