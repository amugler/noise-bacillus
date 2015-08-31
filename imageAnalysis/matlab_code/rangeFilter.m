function [ filter ] = rangeFilter( a, range )
%RANGEFILTER Generate a logical mask within range
%
% a : array to test
% range : two-element array indicating the inclusive range
%         [lower upper]
%
% Mark Kittisopikul, 2010
% UT Southwestern

filter = (a >= range(1)) & (a <= range(2));

end

