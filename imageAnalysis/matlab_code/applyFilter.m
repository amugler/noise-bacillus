function [ outCC, outRP, labelMatrix ] = applyFilter( CC, RP, ranges)
%applyFilter Applies region property filters to connected components and
%regiong properties 

% Adapted from burstfiltering, updateFilters

F = zeros(9,CC.NumObjects);
F(1,:) = rangeFilter([RP.Area],ranges.area);
F(2,:) = rangeFilter([RP.MinIntensity],ranges.min);
F(3,:) = rangeFilter([RP.MeanIntensity],ranges.mean);
F(4,:) = rangeFilter([RP.MaxIntensity],ranges.max);
F(5,:) = rangeFilter([RP.Solidity],ranges.solidity);
F(6,:) = rangeFilter([RP.Extent],ranges.extent);
F(7,:) = rangeFilter([RP.EulerNumber],ranges.euler);
F(8,:) = rangeFilter([RP.Perimeter],ranges.perimeter);
F(9,:) = rangeFilter([RP.Eccentricity],ranges.eccentricity);

filter = all(F,1);
outCC = filtercc(CC,filter);
outRP = RP(filter);
labelMatrix = labelmatrix(outCC);

end

