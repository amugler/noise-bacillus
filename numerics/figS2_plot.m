clear all

dd = '';
fd = '';
fs = 14; fs2 = 12; fs3 = 10; fs4 = 30;
lw = 2; lw2 = 2/2; ms = 20/2; ms2 = 9/2;

load([dd 'native_scan.mat']);
aks_native = aks;
Kstar_native = Kstar;
load([dd 'synex_scan.mat']);
aks_synex = aks;
Kstar_synex = Kstar;

% to convert to dimensionful alpha_k
% native
Gammak = 100;
deltak = .001;
% synex
kk = 10;
lambda = 1e-4;

% transition values
ak1_native = .00072;
ak1_synex = .033;
ak2_native = .0033;
ak2_synex = .7565;

% native:
Kstar_native1 = Kstar_native(1:3,1:29);
[ig,ind] = min(abs(aks_native - ak2_native));
Kstar_native2 = Kstar_native(1,30:ind);
Kstar_native3 = Kstar_native(1,ind+1:end);

% synex:
Kstar_synex1 = Kstar_synex(1:3,1:18);
[ig,ind2] = min(abs(aks_synex - ak2_synex));
Kstar_synex2 = Kstar_synex(1,19:ind2);
Kstar_synex3 = Kstar_synex(1,ind2+1:end);

figure(1); clf

% native
subplot(2,2,1)
h = semilogx(aks_native(1:29)*Gammak*deltak*3600,Kstar_native1(1,:),'b.',...
    aks_native(1:29)*Gammak*deltak*3600,Kstar_native1(2,:),'r.',...
    aks_native(1:29)*Gammak*deltak*3600,Kstar_native1(3,:),'r.',...
    aks_native(30:ind)*Gammak*deltak*3600,Kstar_native2,'r.',...
    aks_native(ind+1:end)*Gammak*deltak*3600,Kstar_native3,'b.',...
    [ak1_native ak1_native]*Gammak*deltak*3600,[0 8],'k--',...
    [ak2_native ak2_native]*Gammak*deltak*3600,[0 8],'k--');
set(h,'linewidth',lw,'markersize',ms)
xlim([.12 10])
ylim([0 8])
xlabel('ComK induction rate, $\alpha_k$ (hr$^{-1}$)',...
    'fontsize',fs,'interpreter','latex')
ylabel('Fixed points, $\bar{n}^*$',...
    'fontsize',fs,'interpreter','latex')
legend(h([1 2]),{'Stable','Unstable'},...
    'fontsize',fs2,'location','se')
set(gca,'fontsize',fs2,'xtick',[.1 1 10])
text(5e-2,8,'A','fontsize',fs4)


% synex
subplot(2,2,2)
h = semilogx(aks_synex(1:18)*kk*lambda*3600,Kstar_synex1(1,:),'b.',...
    aks_synex(1:18)*kk*lambda*3600,Kstar_synex1(2,:),'r.',...
    aks_synex(1:18)*kk*lambda*3600,Kstar_synex1(3,:),'r.',...
    aks_synex(19:ind2)*kk*lambda*3600,Kstar_synex2,'r.',...
    aks_synex(ind2+1:end)*kk*lambda*3600,Kstar_synex3,'b.',...
    [ak1_synex ak1_synex]*kk*lambda*3600,[0 12],'k--',...
    [ak2_synex ak2_synex]*kk*lambda*3600,[0 12],'k--');
set(h,'linewidth',lw,'markersize',ms)
xlim([.04 30])
ylim([0 12])
xlabel('ComK induction rate, $\alpha_k$ (hr$^{-1}$)',...
    'fontsize',fs,'interpreter','latex')
ylabel('Fixed points, $\bar{n}^*$',...
    'fontsize',fs,'interpreter','latex')
legend(h([1 2]),{'Stable','Unstable'},...
    'fontsize',fs2,'location','se')
set(gca,'fontsize',fs2,'xtick',[.1 1 10])
text(9e-3,12,'B','fontsize',fs4)

print(gcf,'-depsc',[fd 'figS2.eps'])
