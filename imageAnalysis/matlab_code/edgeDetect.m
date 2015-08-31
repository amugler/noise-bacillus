function [CC,E] = edgeDetect(I)
% edgeDetect Perform basic edge detection with morphological closing excluding borders
%
% Mark Kittisopikul, 2010
% UT Southwestern

    % edges with laplacian of gaussian method
    E = edge(I,'log',0);
    % close gaps
    E = imclose(E,strel('disk',1));

    % set the outer border as an edge
    E(1:end,1) = 1;
    E(1:end,end) = 1;
    E(1,1:end) = 1;
    E(end,1:end) = 1;

    % get connected component regions
    CC = bwconncomp(~E,4);
end

