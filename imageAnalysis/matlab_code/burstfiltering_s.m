function [ s ] = burstfiltering_s( varargin )
%burstfiltering_s Call burstfiltering but with the output as a struct
%
% Mark Kittisopikul, 2010
% UT Southwestern

[s.ranges,s.cc,s.rp,s.output] = burstfiltering(varargin{:});

end

