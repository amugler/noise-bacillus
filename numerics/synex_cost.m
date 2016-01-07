function [Ssum,kstar,S,n1,n2,c1,c2] = synex_cost(X)

dd = '';
load([dd 'pixel_data.mat'])
load([dd 'synex_scan.mat']);
plotit1 = 0;
plotit2 = 1;
iptg = [0 750 1500 3000 3000]; % nM

% loop over IPTG levels
for i = 1:5
  x = double(synex.px{i})';
  lstr{i} = ['IPTG level ' num2str(i)];

  % experimental CDF with calibration
  if i == 1
    x0 = mode(x);
  end
  x1 = x-x0;
  x1(find(x1 < 0)) = 0;
  n1 = sort(x1/X);
  c1 = (1:length(n1))'/length(n1);
  
  % loop over ak values
  for k = 1:length(aks)
    
    % theoretical CDF
    n2 = (0:N)';
    p2 = sum(pnm(:,:,k),2);
    c2 = cumsum(p2);

    if plotit1
      figure(1); clf
      plot(n1,c1,'b.-',...
           n2,c2,'ro--');
      ylim([-.1 1.1])
      xlabel('ComK molecule number, n')
      ylabel('cumulative probability, c_n')
      title(['IPTG level ' num2str(i) ', ' ...
             'a_k = ' num2str(aks(k)) '; ' ...
             'X = ' num2str(X)])
      drawnow
    end

    % cost function: sum-of-squares per data point
    % 1. find overlap range
    N1 = floor(max(n1));
    N2 = N;
    Nmin = min([N1 N2]);
    % 2. find overlapping points
    dn1 = diff(ceil(n1));
    ind = find(dn1);
    c1_ = [];
    for j = 1:length(ind)
      for z = 1:dn1(ind(j)) % add multiple times if necessary
        c1_ = [c1_; c1(ind(j))];
      end
    end
    c1_ = c1_(1:Nmin);
    c2_ = c2(1+1:Nmin+1);
    % 3. compute sum of squares
    S(i,k) = sum((c1_-c2_).^2)/Nmin;
  end

  % best ak value
  [Smin(i),kstar(i)] = min(S(i,:));
  akstar(i) = aks(kstar(i));
end

if plotit2
  figure(1); clf
  for i = 2:5
    subplot(4,4,i-1)
    x = double(synex.px{i})';
    x1 = x-x0;
    x1(find(x1 < 0)) = 0;
    n1 = sort(x1/X);
    c1 = (1:length(n1))'/length(n1);
    p2 = sum(pnm(:,:,kstar(i)),2);
    c2 = cumsum(p2);
    plot(n1,c1,'b.-',...
         n2,c2,'ro--');
    ylim([-.1 1.1])

    subplot(4,4,i+4-1)
    hold on
    N1 = ceil(max(n1));
    p1 = hist(n1,0:N1);
    p1 = p1/sum(p1);
    N2 = N;
    bar(0:N1,p1);
    plot(0:N2,p2,'r.-');
    xlim([-.5 min([N1 N2])])
    box on
  end

  subplot(4,4,[9 10 13 14])
  semilogx(aks,S,'.-')
  xlabel('a_k')
  ylabel('sum of squares per data point, S')
  legend(lstr,'location','best')

  subplot(4,4,[11 12 15 16])
  plot(iptg,akstar,'ko-')
  xlabel('IPTG (nM)')
  ylabel('a_k*')
  xlim([-100 3100])
  drawnow
end


% sum over IPTG levels
Ssum = sum(Smin);

% output cumulative data
n1 = []; n2 = []; c1 = []; c2 = [];
for i = 1:5
    x = double(synex.px{i})';
    x1 = x-x0;
    x1(find(x1 < 0)) = 0;
    n1{i} = sort(x1/X);
    c1{i} = (1:length(n1{i}))'/length(n1{i});
    n2{i} = (0:N)';
    p2 = sum(pnm(:,:,kstar(i)),2);
    c2{i} = cumsum(p2);
end
