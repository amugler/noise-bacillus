% inflection point method
clear all

dd = '';

% NATIVE

load([dd 'native_scan2.mat']);
aks_native1 = aks;
pns = squeeze(sum(pnm,2));
nvec = (0:N)';

for i = 1:size(pns,2)
  pn = pns(:,i);

  % compute second derivative
  dpn = diff([0;pn]);
  ddpn = diff([0;dpn]);
  
  % smooth out alternations
  ddpn_ = (ddpn + [0;ddpn(1:end-1)])/2;
  
  % find first four inflection points (btw max, "min", and "max")
  ind = find(diff(sign(ddpn_)) ~= 0);

  % unimodal
  if length(ind) < 4
    ind = ind(1:2);

    % center of mode
    i1 = floor((ind(1)+ind(2))/2);

    % width of mode
    i1m = ind(1);
    i1p = ind(2);

    % weight of mode
    pi2_native_inflection1(i) = 1;

    % plot
    clf; hold on; box on
    h1 = fill([nvec(1);nvec;nvec(end)],[0;pn;0],[1 .7 .7]);
    h2 = fill([nvec(i1m);nvec(i1m:i1p);nvec(i1p)],[0;pn(i1m:i1p);0],...
              [1 .3 .3]);
    h3 = plot(nvec(i1),pn(i1),'ko','markersize',8);
    h = plot(nvec,pn,'k.-');

  % "bimodal"
  else
    ind = ind(1:4);
  
    % centers of modes
    i1 = floor((ind(1)+ind(2))/2);
    i2 = floor((ind(3)+ind(4))/2);
    
    % widths of modes
    i1m = ind(1);
    i1p = ind(2);
    i2m = ind(3);
    i2p = ind(4);
    
    % weights of modes
    i0 = ceil((ind(2)+ind(3))/2);
    pi1 = sum(pn(1:i0));
    pi2_native_inflection1(i) = 1-pi1;
    
    % plot
    clf; hold on; box on
    h1 = fill([nvec(1);nvec(1:i0);nvec(i0)],[0;pn(1:i0);0],[.7 .7 1]);
    h2 = fill([nvec(i0);nvec(i0:end);nvec(end)],[0;pn(i0:end);0],...
              [1 .7 .7]);
    h3 = fill([nvec(i1m);nvec(i1m:i1p);nvec(i1p)],[0;pn(i1m:i1p);0],...
              [.3 .3 1]);
    h4 = fill([nvec(i2m);nvec(i2m:i2p);nvec(i2p)],[0;pn(i2m:i2p);0],...
              [1 .3 .3]);
    h5 = plot(nvec(i1),pn(i1),'ko','markersize',8);
    h5 = plot(nvec(i2),pn(i2),'ko','markersize',8);
    h = plot(nvec,pn,'k.-');
  end
  pmax = max(pn);
  ylim([-.1*pmax pmax*1.1])
  drawnow; pause(.02)
end

load([dd 'native_scan.mat']);
aks_native2 = aks;
pns = squeeze(sum(pnm,2));
nvec = (0:N)';

for i = 1:size(pns,2)
  pn = pns(:,i);

  % compute second derivative
  dpn = diff([0;pn]);
  ddpn = diff([0;dpn]);
  
  % smooth out alternations
  ddpn_ = (ddpn + [0;ddpn(1:end-1)])/2;
  
  % find first four inflection points (btw max, "min", and "max")
  ind = find(diff(sign(ddpn_)) ~= 0);

  % unimodal
  if length(ind) < 4
    ind = ind(1:2);

    % center of mode
    i1 = floor((ind(1)+ind(2))/2);

    % width of mode
    i1m = ind(1);
    i1p = ind(2);

    % weight of mode
    pi2_native_inflection2(i) = 1;

    % plot
    clf; hold on; box on
    h1 = fill([nvec(1);nvec;nvec(end)],[0;pn;0],[1 .7 .7]);
    h2 = fill([nvec(i1m);nvec(i1m:i1p);nvec(i1p)],[0;pn(i1m:i1p);0],...
              [1 .3 .3]);
    h3 = plot(nvec(i1),pn(i1),'ko','markersize',8);
    h = plot(nvec,pn,'k.-');

  % "bimodal"
  else
    ind = ind(1:4);
  
    % centers of modes
    i1 = floor((ind(1)+ind(2))/2);
    i2 = floor((ind(3)+ind(4))/2);
    
    % widths of modes
    i1m = ind(1);
    i1p = ind(2);
    i2m = ind(3);
    i2p = ind(4);
    
    % weights of modes
    i0 = ceil((ind(2)+ind(3))/2);
    pi1 = sum(pn(1:i0));
    pi2_native_inflection2(i) = 1-pi1;
    
    % plot
    clf; hold on; box on
    h1 = fill([nvec(1);nvec(1:i0);nvec(i0)],[0;pn(1:i0);0],[.7 .7 1]);
    h2 = fill([nvec(i0);nvec(i0:end);nvec(end)],[0;pn(i0:end);0],...
              [1 .7 .7]);
    h3 = fill([nvec(i1m);nvec(i1m:i1p);nvec(i1p)],[0;pn(i1m:i1p);0],...
              [.3 .3 1]);
    h4 = fill([nvec(i2m);nvec(i2m:i2p);nvec(i2p)],[0;pn(i2m:i2p);0],...
              [1 .3 .3]);
    h5 = plot(nvec(i1),pn(i1),'ko','markersize',8);
    h5 = plot(nvec(i2),pn(i2),'ko','markersize',8);
    h = plot(nvec,pn,'k.-');
  end
  pmax = max(pn);
  ylim([-.1*pmax pmax*1.1])
  drawnow; pause(.02)
end

aks_native = [aks_native1 aks_native2];
pi2_native_inflection = [pi2_native_inflection1 pi2_native_inflection2];

% SYNEX

load([dd 'synex_scan2.mat']);
aks_synex1 = aks;
pns = squeeze(sum(pnm,2));
nvec = (0:N)';

for i = 1:size(pns,2)
  pn = pns(:,i);

  % compute second derivative
  dpn = diff([0;pn]);
  ddpn = diff([0;dpn]);

  % smooth out alternations
  ddpn_ = (ddpn + [0;ddpn(1:end-1)])/2;

  % find first four inflection points (btw max, "min", and "max")
  ind = find(diff(sign(ddpn_)) ~= 0);

  % unimodal
  if length(ind) < 4
    ind = ind(1:2);

    % center of mode
    i1 = floor((ind(1)+ind(2))/2);

    % width of mode
    i1m = ind(1);
    i1p = ind(2);

    % weight of mode
    pi2_synex_inflection1(i) = 1;

    % plot
    clf; hold on; box on
    h1 = fill([nvec(1);nvec;nvec(end)],[0;pn;0],[1 .7 .7]);
    h2 = fill([nvec(i1m);nvec(i1m:i1p);nvec(i1p)],[0;pn(i1m:i1p);0],...
              [1 .3 .3]);
    h3 = plot(nvec(i1),pn(i1),'ko','markersize',8);
    h = plot(nvec,pn,'k.-');

  % "bimodal"
  else
    ind = ind(1:4);

    % centers of modes
    i1 = floor((ind(1)+ind(2))/2);
    i2 = floor((ind(3)+ind(4))/2);

    % widths of modes
    i1m = ind(1);
    i1p = ind(2);
    i2m = ind(3);
    i2p = ind(4);

    % weights of modes
    i0 = ceil((ind(2)+ind(3))/2);
    pi1 = sum(pn(1:i0));
    pi2_synex_inflection1(i) = 1-pi1;

    % plot
    clf; hold on; box on
    h1 = fill([nvec(1);nvec(1:i0);nvec(i0)],[0;pn(1:i0);0],[.7 .7 ...
                        1]);
    h2 = fill([nvec(i0);nvec(i0:end);nvec(end)],[0;pn(i0:end);0],...
              [1 .7 .7]);
    h3 = fill([nvec(i1m);nvec(i1m:i1p);nvec(i1p)],[0;pn(i1m:i1p);0],...
              [.3 .3 1]);
    h4 = fill([nvec(i2m);nvec(i2m:i2p);nvec(i2p)],[0;pn(i2m:i2p);0],...
              [1 .3 .3]);
    h5 = plot(nvec(i1),pn(i1),'ko','markersize',8);
    h5 = plot(nvec(i2),pn(i2),'ko','markersize',8);
    h = plot(nvec,pn,'k.-');
  end
  pmax = max(pn);
  ylim([-.1*pmax pmax*1.1])
  drawnow; pause(.02)
end

load([dd 'synex_scan.mat']);
aks_synex2 = aks;
pns = squeeze(sum(pnm,2));
nvec = (0:N)';

for i = 1:size(pns,2)
  pn = pns(:,i);

  % compute second derivative
  dpn = diff([0;pn]);
  ddpn = diff([0;dpn]);

  % smooth out alternations
  ddpn_ = (ddpn + [0;ddpn(1:end-1)])/2;

  % find first four inflection points (btw max, "min", and "max")
  ind = find(diff(sign(ddpn_)) ~= 0);

  % unimodal
  if length(ind) < 4
    ind = ind(1:2);

    % center of mode
    i1 = floor((ind(1)+ind(2))/2);

    % width of mode
    i1m = ind(1);
    i1p = ind(2);

    % weight of mode
    pi2_synex_inflection2(i) = 1;

    % plot
    clf; hold on; box on
    h1 = fill([nvec(1);nvec;nvec(end)],[0;pn;0],[1 .7 .7]);
    h2 = fill([nvec(i1m);nvec(i1m:i1p);nvec(i1p)],[0;pn(i1m:i1p);0],...
              [1 .3 .3]);
    h3 = plot(nvec(i1),pn(i1),'ko','markersize',8);
    h = plot(nvec,pn,'k.-');

  % "bimodal"
  else
    ind = ind(1:4);

    % centers of modes
    i1 = floor((ind(1)+ind(2))/2);
    i2 = floor((ind(3)+ind(4))/2);

    % widths of modes
    i1m = ind(1);
    i1p = ind(2);
    i2m = ind(3);
    i2p = ind(4);

    % weights of modes
    i0 = ceil((ind(2)+ind(3))/2);
    pi1 = sum(pn(1:i0));
    pi2_synex_inflection2(i) = 1-pi1;

    % plot
    clf; hold on; box on
    h1 = fill([nvec(1);nvec(1:i0);nvec(i0)],[0;pn(1:i0);0],[.7 .7 ...
                        1]);
    h2 = fill([nvec(i0);nvec(i0:end);nvec(end)],[0;pn(i0:end);0],...
              [1 .7 .7]);
    h3 = fill([nvec(i1m);nvec(i1m:i1p);nvec(i1p)],[0;pn(i1m:i1p);0],...
              [.3 .3 1]);
    h4 = fill([nvec(i2m);nvec(i2m:i2p);nvec(i2p)],[0;pn(i2m:i2p);0],...
              [1 .3 .3]);
    h5 = plot(nvec(i1),pn(i1),'ko','markersize',8);
    h5 = plot(nvec(i2),pn(i2),'ko','markersize',8);
    h = plot(nvec,pn,'k.-');
  end
  pmax = max(pn);
  ylim([-.1*pmax pmax*1.1])
  drawnow; pause(.02)
end

aks_synex = [aks_synex1 aks_synex2];
pi2_synex_inflection = [pi2_synex_inflection1 pi2_synex_inflection2];

% SAVE
save([dd 'pi2_inflection_full.mat'],'pi2_native_inflection',...
     'pi2_synex_inflection','aks_native','aks_synex')