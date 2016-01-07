clear all

dd = '../dat/';
fd = '../fig/';

% synex
iptg = [0 750 1500 3000 30e3]; % nM
load([dd 'pixel_data.mat'])
for i = 1:4
  xs{i} = double(synex.px{i})';
end
load([dd 'pixels_only_30um_2013_07_31.mat'])
xs{5} = double(synex.induced30um);

% native
iptg = [0 750 1500 3000 30e3 100e3]; % nM
load([dd 'pixel_data.mat'])
for i = 1:4
  xn{i} = double(native.px{i})';
end
load([dd 'pixels_only_30um_2013_07_31.mat'])
xn{5} = double(native.induced30um);
load([dd 'M1240.mat'])
xn{6} = allpx;

% histograms
M = 0:10:1000;
for i = 1:4
  [cs{i},fs{i}] = hist(xs{i},M);
end
M = 0:20:2000;
[cs{5},fs{5}] = hist(xs{5},M);

M = 0:10:1000;
for i = 1:4
  [cn{i},fn{i}] = hist(xn{i},M);
end
M = 0:20:2000;
[cn{5},fn{5}] = hist(xn{5},M);
M = 0:100:10000;
[cn{6},fn{6}] = hist(xn{6},M);


% plot
figure(1); clf
lw = .1; 
for i = 2:4
  subplot(5,5,i+5-1)
  hs(i-1) = bar(fs{i},cs{i}/1e3,'r');
  xlim([180 500])
  ylim([0 max(cs{i}/1e3)*1.1])
  set(gca,'xtick',[300 500])
end
subplot(5,5,5+5-1)
hs(5-1) = bar(fs{5},cs{5}/1e3,'r');
xlim([180 1000])
ylim([0 max(cs{5}/1e3)*1.1])
set(gca,'xtick',[500 1000])
for i = 2:4
  subplot(5,5,i-1)
  hn(i-1) = bar(fn{i},cn{i}/1e3,'b');
  xlim([180 500])
  ylim([0 max(cn{i}/1e3)*1.1])
  set(gca,'xtick',[300 500])
end
subplot(5,5,5-1)
hn(5-1) = bar(fn{5},cn{5}/1e3,'b');
xlim([180 1000])
ylim([0 max(cn{5}/1e3)*1.1])
set(gca,'xtick',[500 1000])
subplot(5,5,6-1)
hn(6-1) = bar(fn{6},cn{6}/1e3,'b');
xlim([180 7000])
ylim([0 max(cn{6}/1e3)*1.1])
set(gca,'xtick',[2000 6000])

subplot(5,5,1); ylabel('count / 10^3')
subplot(5,5,6); ylabel('count / 10^3')
for i = 5:9
  subplot(5,5,i)
  xlabel('fluorescence (AU)')
end
set(hs,'linewidth',lw,'edgecolor','r')
set(hn,'linewidth',lw,'edgecolor','b')

print(gcf,'-depsc',[fd 'fig3b_v2.eps'])