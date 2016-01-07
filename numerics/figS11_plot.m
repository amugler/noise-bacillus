clear all

dd = '';
fd = '';

% synex
iptg = [0 0 3e3]; % nM
load([dd 'pixel_data.mat'])
xs{1} = double(synex.px{1})';
xs{3} = double(synex.px{5})';
load([dd 'pixels_only_30um_2013_07_31.mat'])
xs{2} = double(synex.not_induced);

% native
iptg = [0 0 30e3]; % nM
load([dd 'pixel_data.mat'])
xn{1} = double(native.px{1})';
load([dd 'pixels_only_30um_2013_07_31.mat'])
xn{2} = double(native.not_induced);
xn{3} = double(native.induced30um_later);

% histograms
M = 0:10:1000;
[cs{1},fs{1}] = hist(xs{1},M);
M = 0:10:1000;
[cs{2},fs{2}] = hist(xs{2},M);
M = 0:20:2000;
[cs{3},fs{3}] = hist(xs{3},M);

M = 0:10:1000;
[cn{1},fn{1}] = hist(xn{1},M);
M = 0:10:1000;
[cn{2},fn{2}] = hist(xn{2},M);
M = 0:20:2000;
[cn{3},fn{3}] = hist(xn{3},M);

% plot
figure(1); clf
lw = .1; 

subplot(10,5,[16 21])
hs(1) = bar(fs{1},cs{1}/1e3,'r');
xlim([180 500])
ylim([0 max(cs{1}/1e3)*1.1])
set(gca,'xtick',[300 500])
title('[IPTG] = 0 $\mu$M','interpreter','latex')

subplot(10,5,[17 22])
hs(2) = bar(fs{2},cs{2}/1e3,'r');
xlim([180 500])
ylim([0 max(cs{2}/1e3)*1.1])
set(gca,'xtick',[300 500])
title('[IPTG] = 0 $\mu$M','interpreter','latex')

subplot(10,5,[18 23])
hs(3) = bar(fs{3},cs{3}/1e3,'r');
xlim([180 750])
ylim([0 max(cs{3}/1e3)*1.1])
set(gca,'xtick',[250 500 750])
title('[IPTG] = 3 $\mu$M','interpreter','latex')

subplot(10,5,[1 6])
hn(1) = bar(fn{1},cn{1}/1e3,'b');
xlim([180 500])
ylim([0 max(cn{1}/1e3)*1.1])
set(gca,'xtick',[300 500])
title('[IPTG] = 0 $\mu$M','interpreter','latex')

subplot(10,5,[2 7])
hn(2) = bar(fn{2},cn{2}/1e3,'b');
xlim([180 500])
ylim([0 max(cn{2}/1e3)*1.1])
set(gca,'xtick',[300 500])
title('[IPTG] = 0 $\mu$M','interpreter','latex')

subplot(10,5,[3 8])
hn(3) = bar(fn{3},cn{3}/1e3,'b');
xlim([180 1200])
ylim([0 max(cn{3}/1e3)*1.1])
set(gca,'xtick',[600 1200])
title('[IPTG] = 30 $\mu$M','interpreter','latex')

set(hs,'linewidth',lw,'edgecolor','r')
set(hn,'linewidth',lw,'edgecolor','b')

print(gcf,'-depsc',[fd 'figS11.eps'])