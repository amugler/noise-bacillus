function rgb = makergb(varargin)
% Construct a rgb matrix where each color channel is normalized
%
% INPUT
% r - red channel with dimensions Y by X
% g - green channel with dimensions Y by X
% b - blue channel with dimensions Y by X
% ... additional channels
% normalizeFlag - true if each channel should be normalized (default)
%               - false otherwise
%
% OUTPUT
% rgb a uint8 Y by X by N where N is at least 3 and corresponds with the
% number of channels given in the input. Empty matrices are expanded to
% zeros.
%
% See also mat2gray
%
% Mark Kittisopikul

    % size of the red channel
    siz = size(varargin{1});
    % number of channels
    N = nargin;

    % detect normalization flag, which should be a scalar as the last input
    normalize = true;
    if(isscalar(varargin{N}))
        N = N-1;
        normalize = logical(varargin{N});
    end

    % ensure that the number of channels is at least 3
    [varargin{N+1:3}] = deal(zeros(siz));
    N = 3;

    % convert any empty channels to zeros
    empty = cellfun('isempty',varargin(1:N));
    [varargin{empty}] = deal(zeros(siz));

    % normalize each channel
    if(normalize)
        varargin = cellfun(@mat2gray,varargin(1:N),'UniformOutput',false);
    end

    % assemble final matrix and convert to rgb
    rgb = cat(3,varargin{1:N});
    rgb = im2uint8(rgb);

end