clear all
dd = '';
fd = '';
fs = 25; fs2 = 15; lw = 3; lw2 = 2; ms = 20; ms2 = 9;

load([dd 'native_scan.mat']);
aks_native = aks;
load([dd 'synex_scan.mat']);
aks_synex = aks;
load([dd 'pi2_inflection.mat']);
load([dd 'pi2_dkl.mat']);

ak2_native = .0033;
ak2_synex = .7565;
akbar_native = aks_native/ak2_native;
akbar_synex = aks_synex/ak2_synex;

figure(1); clf
h = semilogx(akbar_native,pi2_native_inflection,'b-',...
    akbar_native,pi2_native_dkl,'b--',...
    akbar_synex,pi2_synex_inflection,'r-',...
    akbar_synex,pi2_synex_dkl,'r--',...
    [1 1],[-.05 1.05],'k--');
ylim([-.05 1.05])
xlim([2e-2 4e1])
set(h,'linewidth',lw)
legend(h,{'Native: inflection points',...
    'Native: two Poissons',...
    'SynEx: inflection points',...
    'SynEx: two Poissons',...
    'Deterministic transition'},...
    'location','nw','fontsize',fs2)
xlabel('ComK induction rate, $a_k$ (normalized)',...
    'fontsize',fs,'interpreter','latex')
ylabel('Fraction in responsive state, $f$',...
    'fontsize',fs,'interpreter','latex')
set(gca,'fontsize',fs2)
box on
print(gcf,'-depsc',[fd 'figS4.eps'])