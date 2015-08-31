function [ newcc, oldccidx ] = filtercc( cc, logicfilter )
%filtercc Perform logical and index filtering on connected components structure
%
% Mark Kittisopikul, 2010
% UT Southwestern

newcc = cc;
newcc.PixelIdxList = cc.PixelIdxList(logicfilter);
newcc.NumObjects = length(newcc.PixelIdxList);
idx = 1:cc.NumObjects;
oldccidx = idx(logicfilter);

end

