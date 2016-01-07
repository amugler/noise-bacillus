% synex: intermediate induction
% computes error bar on a_k for experimental data
clear all

dd = '';
load([dd 'synex_scan.mat'])
load([dd 'synex_fit_cdf2.mat']);

f = .25;
for i = 1:2
  [Smin,kstar] = min(S(i,:));
  Scross(i,1) = Smin*(1+f);
  dS = S(i,:)-Scross(i);
  if dS(1) < 0
    kcross{i} =[1 find(diff(sign(dS)))];
  elseif dS(end) < 0
    kcross{i} = [find(diff(sign(dS))) length(dS)];
  else
    kcross{i} = find(diff(sign(dS)));
  end
end

figure(1); clf;
loglog(aks,S,'.-',...
       [min(aks) max(aks)],[Scross(1) Scross(1)],'-',...
       [min(aks) max(aks)],[Scross(2) Scross(2)],'-')

save([dd 'error_synex2.mat'],'kcross')
