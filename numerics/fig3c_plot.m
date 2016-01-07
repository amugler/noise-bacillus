clear all

dd = '';
fd = '';
fs = 14; fs2 = 12; fs3 = 10; lw = 2; lw2 = 2/2; ms = 20/2; ms2 = 9/2;

load([dd 'pi2_inflection_full.mat']);
load([dd 'pi2_dkl_full.mat']);

% to convert to dimensionful alpha_k
% native
Gammak = 100;
deltak = .001;
% synex
kk = 10;
lambda = 1e-4;

% deterministic values
ak1_native = .00072;
ak1_synex = .033;
ak2_native = .0033;
ak2_synex = .7565;

% inflection smoother at high/low stress
% dkl smoother at interediate stress
k1n = 40;
k2n = min(intersect(find(pi2_native_inflection == 1),...
    k1n+1:length(pi2_native_inflection)));
aks_native = aks_native([1:k2n k2n:end]);
pi2_native_theory = [pi2_native_inflection(1:k1n) ...
    pi2_native_dkl(k1n+1:k2n) ...
    ones(1,length(aks_native)-k2n)];
k1s = 40;
k2s = min(intersect(find(pi2_synex_inflection == 1),...
    k1s+1:length(pi2_synex_inflection)));
aks_synex = aks_synex([1:k2s k2s:end]);
pi2_synex_theory = [pi2_synex_inflection(1:k1s) ...
    pi2_synex_dkl(k1s+1:k2s) ...
    ones(1,length(aks_synex)-k2s)];

% native: first 24 points are unimodal and (mistakenly)
%   identified as pi2 = 1 instead of 0 due to symmetry
pi2_native_theory(1:24) = 0;

% use smoothing at low-stress stitch point
T = 5; % number of smoothing iterations
di = 2; % shift by di each iteration
N = 5; % smoothing distance
M = 10; % smoothing region
cn = round(k1n - di*(T-1)/2);
cs = round(k1n - di*(T-1)/2);
for t = 1:T
    for i = cn-M:cn+M
        pi2_native_theory(i) = mean(pi2_native_theory(i-N:i+N));
    end
    for i = cs-M:cs+M
        pi2_synex_theory(i) = mean(pi2_synex_theory(i-N:i+N));
    end
    cn = cn + di;
    cs = cs + di;
end

% expansion factors
fmin = 1e-2;
ak0_native = aks_native(min(find(pi2_native_theory > fmin)));
ak3_native = aks_native(min(find(pi2_native_theory == 1)));
ak0_synex = aks_synex(min(find(pi2_synex_theory > fmin)));
ak3_synex = aks_synex(min(find(pi2_synex_theory == 1)));

alphak0_native = ak0_native*Gammak*deltak*3600;
alphak1_native = ak1_native*Gammak*deltak*3600;
alphak2_native = ak2_native*Gammak*deltak*3600;
alphak3_native = ak3_native*Gammak*deltak*3600;
alphak0_synex = ak0_synex*kk*lambda*3600;
alphak1_synex = ak1_synex*kk*lambda*3600;
alphak2_synex = ak2_synex*kk*lambda*3600;
alphak3_synex = ak3_synex*kk*lambda*3600;

exp1_native = ak1_native/ak0_native
exp2_native = ak3_native/ak2_native
exp1_synex = ak1_synex/ak0_synex
exp2_synex = ak3_synex/ak2_synex

exp_total_native = (ak3_native/ak0_native)/(ak2_native/ak1_native)
exp_total_synex = (ak3_synex/ak0_synex)/(ak2_synex/ak1_synex)

figure(1); clf

% native
subplot(2,2,1)
h = semilogx(aks_native*Gammak*deltak*3600,pi2_native_theory,'b-',...
    [ak1_native ak1_native]*Gammak*deltak*3600,[-.05 1.05],'k--',...
    [ak2_native ak2_native]*Gammak*deltak*3600,[-.05 1.05],'k--');
set(h,'linewidth',lw)
xlim([min(aks_native)*Gammak*deltak*3600 ...
    max(aks_native)*Gammak*deltak*3600])
ylim([-.05 1.05])
xlabel('ComK induction rate, $\alpha_k$ (hr$^{-1}$)',...
    'fontsize',fs,'interpreter','latex')
ylabel('Frac.\ in responsive state, $f$',...
    'fontsize',fs,'interpreter','latex')
text(7e-5*Gammak*deltak*3600,.65,'excitable',...
    'fontsize',fs3,'horizontalalignment','left')
text(sqrt(ak1_native*ak2_native)*Gammak*deltak*3600,.65,'osc.',...
    'fontsize',fs3,'horizontalalignment','center')
text(4.5e-3*Gammak*deltak*3600,.65,'mono-',...
    'fontsize',fs3,'horizontalalignment','left')
text(4.5e-3*Gammak*deltak*3600,.55,'stable',...
    'fontsize',fs3,'horizontalalignment','left')
legend(h([3 1]),{'deterministic','stochastic'},...
    'fontsize',fs2,'location','nw')
set(gca,'fontsize',fs2,'xtick',[.001 .01 .1 1 10])


% synex
subplot(2,2,3)
h = semilogx(aks_synex*kk*lambda*3600,pi2_synex_theory,'r-',...
    [ak1_synex ak1_synex]*kk*lambda*3600,[-.05 1.05],'k--',...
    [ak2_synex ak2_synex]*kk*lambda*3600,[-.05 1.05],'k--','linewidth',lw);
set(h,'linewidth',lw)
xlim([min(aks_synex)*kk*lambda*3600 ...
    max(aks_synex)*kk*lambda*3600])
ylim([-.05 1.05])
xlabel('ComK induction rate, $\alpha_k$ (hr$^{-1}$)',...
    'fontsize',fs,'interpreter','latex')
ylabel('Frac.\ in responsive state, $f$',...
    'fontsize',fs,'interpreter','latex')
text(3e-3*kk*lambda*3600,.65,'excitable',...
    'fontsize',fs3,'horizontalalignment','left')
text(sqrt(ak1_synex*ak2_synex)*kk*lambda*3600,.65,'oscillatory',...
    'fontsize',fs3,'horizontalalignment','center')
text(1*kk*lambda*3600,.65,'m.s.',...
    'fontsize',fs3,'horizontalalignment','left')
legend(h([3 1]),{'deterministic','stochastic'},...
    'fontsize',fs2,'location','nw')
set(gca,'fontsize',fs2,'xtick',[.001 .01 .1 1 10])

print(gcf,'-depsc',[fd 'fig3c.eps'])
