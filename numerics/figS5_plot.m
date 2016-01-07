clear all
dd = '';
fd = '';
fs = 12; fs2 = 10; fs3 = 30; lw = 1.5; lw2 = 1; ms = 10; ms2 = 5;

load([dd 'native_scan.mat']);
aks_A = aks;
Gammak = 100;
deltak = .001;
load([dd 'synex_scan.mat']);
aks_B = aks;
kk = 10;
lambda = 1e-4;
load([dd 'native_fit_cdf.mat']);
n1_1 = n1{3};
n2_1 = n2{3};
c1_1 = c1{3};
c2_1 = c2{3};
S_1 = S(3,:);
n1_2 = n1{4};
n2_2 = n2{4};
c1_2 = c1{4};
c2_2 = c2{4};
S_2 = S(4,:);
load([dd 'native_fit_cdf3.mat']);
n1_3 = n1{1};
n2_3 = n2{1};
c1_3 = c1{1};
c2_3 = c2{1};
S_3 = S(1,:);
load([dd 'synex_fit_cdf.mat']);
n1_4 = n1{2};
n2_4 = n2{2};
c1_4 = c1{2};
c2_4 = c2{2};
S_4 = S(2,:);
n1_5 = n1{3};
n2_5 = n2{3};
c1_5 = c1{3};
c2_5 = c2{3};
S_5 = S(3,:);
load([dd 'synex_fit_cdf2.mat']);
n1_6 = n1{2};
n2_6 = n2{2};
c1_6 = c1{2};
c2_6 = c2{2};
S_6 = S(2,:);

figure(1); clf
subplot(2,2,1)
h = plot(n1_1,c1_1,'b-',n2_1,c2_1,'r.-',...
    'linewidth',lw,'markersize',ms);
xlim([0 80])
ylim([0 1.1])
ylabel('Cumulative probability','fontsize',fs)
xlabel('ComK molecule #','fontsize',fs)
title('[IPTG] = 1.5 \muM')
legend(h,{'Experiment','Theory'},'fontsize',fs2,'location','se')
set(gca,'fontsize',fs2)
text(-15,1.1,'A','fontsize',fs3)

subplot(2,2,2)
h = plot(n1_2,c1_2,'b-',n2_2,c2_2,'r.-',...
    'linewidth',lw,'markersize',ms);
xlim([0 80])
ylim([0 1.1])
ylabel('Cumulative probability','fontsize',fs)
xlabel('ComK molecule #','fontsize',fs)
title('[IPTG] = 3 \muM')
legend(h,{'Experiment','Theory'},'fontsize',fs2,'location','se')
set(gca,'fontsize',fs2)
text(-15,1.1,'B','fontsize',fs3)

subplot(2,2,3)
h = plot(n1_3,c1_3,'b-',n2_3,c2_3,'r.-',...
    'linewidth',lw,'markersize',ms);
xlim([0 80])
ylim([0 1.1])
ylabel('Cumulative probability','fontsize',fs)
xlabel('ComK molecule #','fontsize',fs)
title('[IPTG] = 100 \muM')
legend(h,{'Experiment','Theory'},'fontsize',fs2,'location','se')
set(gca,'fontsize',fs2)
text(-15,1.1,'C','fontsize',fs3)

subplot(2,2,4)
h = semilogx(aks_A*Gammak*deltak*3600,S_1,'k.-',...
    aks_A*Gammak*deltak*3600,S_2,'g.-',...
    aks_A*Gammak*deltak*3600,S_3,'m.-',...
    'linewidth',lw,'markersize',ms);
xlim([min(aks_A)*Gammak*deltak*3600 ...
    max(aks_A)*Gammak*deltak*3600])
xlabel('ComK induction rate, $\alpha_k$ (hr$^{-1}$)',...
    'fontsize',fs,'interpreter','latex')
ylabel('Sum of squares, $S$',...
    'fontsize',fs,'interpreter','latex')
legend(h,{'[IPTG] = 1.5 \muM','[IPTG] = 3 \muM','[IPTG] = 100 \muM'},...
    'fontsize',fs2,'location','nw')
set(gca,'fontsize',fs2,'xtick',[.1 1 10])
text(8e-3,.6,'D','fontsize',fs3)

print(gcf,'-depsc',[fd 'figS5-1.eps'])

figure(2); clf
subplot(2,2,1)
h = plot(n1_4,c1_4,'b-',n2_4,c2_4,'r.-',...
    'linewidth',lw,'markersize',ms);
xlim([0 50])
ylim([0 1.1])
ylabel('Cumulative probability','fontsize',fs)
xlabel('ComK molecule #','fontsize',fs)
title('[IPTG] = 0.75 \muM')
legend(h,{'Experiment','Theory'},'fontsize',fs2,'location','se')
set(gca,'fontsize',fs2)
text(-10,1.1,'E','fontsize',fs3)

subplot(2,2,2)
h = plot(n1_5,c1_5,'b-',n2_5,c2_5,'r.-',...
    'linewidth',lw,'markersize',ms);
xlim([0 50])
ylim([0 1.1])
ylabel('Cumulative probability','fontsize',fs)
xlabel('ComK molecule #','fontsize',fs)
title('[IPTG] = 1.5 \muM')
legend(h,{'Experiment','Theory'},'fontsize',fs2,'location','se')
set(gca,'fontsize',fs2)
text(-10,1.1,'F','fontsize',fs3)

subplot(2,2,3)
h = plot(n1_6,c1_6,'b-',n2_6,c2_6,'r.-',...
    'linewidth',lw,'markersize',ms);
xlim([0 50])
ylim([0 1.1])
ylabel('Cumulative probability','fontsize',fs)
xlabel('ComK molecule #','fontsize',fs)
title('[IPTG] = 3 \muM')
legend(h,{'Experiment','Theory'},'fontsize',fs2,'location','se')
set(gca,'fontsize',fs2)
text(-10,1.1,'G','fontsize',fs3)

subplot(2,2,4)
h = semilogx(aks_B*kk*lambda*3600,S_4,'k.-',...
    aks_B*kk*lambda*3600,S_5,'g.-',...
    aks_B*kk*lambda*3600,S_6,'m.-',...
    'linewidth',lw,'markersize',ms);
xlim([min(aks_B)*kk*lambda*3600 ...
    max(aks_B)*kk*lambda*3600])
xlabel('ComK induction rate, $\alpha_k$ (hr$^{-1}$)',...
    'fontsize',fs,'interpreter','latex')
ylabel('Sum of squares, $S$',...
    'fontsize',fs,'interpreter','latex')
legend(h,{'[IPTG] = 0.75 \muM','[IPTG] = 1.5 \muM','[IPTG] = 3 \muM'},...
    'fontsize',fs2,'location','nw')
set(gca,'fontsize',fs2,'xtick',[.1 1 10])
text(8e-3,.4,'H','fontsize',fs3)

print(gcf,'-depsc',[fd 'figS5-2.eps'])
