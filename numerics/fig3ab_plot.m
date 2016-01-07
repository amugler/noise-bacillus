clear all

dd = '';
fd = '';
ms = 30/3; lw = 5/3; lw2 = 2/3; fs = 36/3; fs2 = 30/3; fs3 = 20/3;
c1 = [.5 .5 1]; c2 = [1 .5 .5];

% Fig 3A: time series
figure(1); clf

% native
subplot(3,3,1)
load([dd 'native_excitable.mat'])
plot(t/3600,K,'b-','linewidth',lw)
xlim([0 40])
ylim([-3 50])
xlabel('Time (hours)','fontsize',fs)
ylabel('ComK molecule #','fontsize',fs)
title('Excitable','fontsize',fs)
set(gca,'fontsize',fs2)

subplot(3,3,2)
load([dd 'native_oscillatory.mat'])
plot(t/3600,K,'b-','linewidth',lw)
xlim([0 40])
ylim([2 8])
xlabel('Time (hours)','fontsize',fs)
ylabel('ComK molecule #','fontsize',fs)
title('Oscillatory','fontsize',fs)
set(gca,'fontsize',fs2)

subplot(3,3,3)
load([dd 'native_monostable.mat'])
plot(t/3600,K,'b-','linewidth',lw)
xlim([0 40])
ylim([0 60])
xlabel('Time (hours)','fontsize',fs)
ylabel('ComK molecule #','fontsize',fs)
title('Mono-stable','fontsize',fs)
set(gca,'fontsize',fs2)

% synex
subplot(3,3,4)
load([dd 'synex_excitable.mat'])
plot(t/3600,K,'r-','linewidth',lw)
xlim([0 40])
ylim([-3 25])
xlabel('Time (hours)','fontsize',fs)
ylabel('ComK molecule #','fontsize',fs)
set(gca,'fontsize',fs2)

subplot(3,3,5)
load([dd 'synex_oscillatory.mat'])
plot(t/3600,K,'r-','linewidth',lw)
xlim([0 40])
ylim([-3 25])
xlabel('Time (hours)','fontsize',fs)
ylabel('ComK molecule #','fontsize',fs)
set(gca,'fontsize',fs2)

subplot(3,3,6)
load([dd 'synex_monostable.mat'])
plot(t/3600,K,'r-','linewidth',lw)
xlim([0 40])
ylim([-3 25])
xlabel('Time (hours)','fontsize',fs)
ylabel('ComK molecule #','fontsize',fs)
set(gca,'fontsize',fs2)

print(gcf,'-depsc',[fd 'fig3a.eps'])


% Fig 3B: distributions
figure(2); clf

% native
subplot(3,3,1)
load([dd 'native_excitable.mat'])
plot(nvec,pn,'b.-','linewidth',lw,'markersize',ms)
xlim([0 80])
xlabel('ComK molecule #','fontsize',fs)
ylabel('Probability','fontsize',fs)
set(gca,'fontsize',fs2)

subplot(3,3,2); hold on
load([dd 'native_oscillatory.mat'])
fill([nvec(i0);nvec(i0:end);nvec(end)],[0;pn(i0:end);0],c1);
plot(nvec,pn,'b.-','linewidth',lw,'markersize',ms)
xlim([0 80])
ylim([0 .17])
xlabel('ComK molecule #','fontsize',fs)
ylabel('Probability','fontsize',fs)
set(gca,'fontsize',fs2)
box on;

subplot(3,3,3); hold on
load([dd 'native_monostable.mat'])
fill([nvec(1);nvec;nvec(end)],[0;pn;0],c1);
plot(nvec,pn,'b.-','linewidth',lw,'markersize',ms)
xlim([0 80])
xlabel('ComK molecule #','fontsize',fs)
ylabel('Probability','fontsize',fs)
set(gca,'fontsize',fs2)
box on

% synex
subplot(3,3,4)
load([dd 'synex_excitable.mat'])
plot(nvec,pn,'r.-','linewidth',lw,'markersize',ms)
xlim([0 80])
xlabel('ComK molecule #','fontsize',fs)
ylabel('Probability','fontsize',fs)
set(gca,'fontsize',fs2)

subplot(3,3,5); hold on
load([dd 'synex_oscillatory.mat'])
fill([nvec(i0);nvec(i0:end);nvec(end)],[0;pn(i0:end);0],c2);
plot(nvec,pn,'r.-','linewidth',lw,'markersize',ms)
xlim([0 80])
ylim([0 .17])
xlabel('ComK molecule #','fontsize',fs)
ylabel('Probability','fontsize',fs)
set(gca,'fontsize',fs2)
box on

subplot(3,3,6); hold on
load([dd 'synex_monostable.mat'])
fill([nvec(1);nvec;nvec(end)],[0;pn;0],c2);
plot(nvec,pn,'r.-','linewidth',lw,'markersize',ms)
xlim([0 80])
xlabel('ComK molecule #','fontsize',fs)
ylabel('Probability','fontsize',fs)
set(gca,'fontsize',fs2)
box on

print(gcf,'-depsc',[fd 'fig3b.eps'])
