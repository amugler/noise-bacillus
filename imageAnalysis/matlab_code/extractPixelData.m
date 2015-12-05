function [ native,synex ] = extractPixelData( rangesFile, imgBaseDir, manual )
%extractPixelData Extract pixel data from images for histograms

if(nargout < 1)
    path = fileparts(mfilename('fullpath'));
    rangesFile = [ path filesep 'regionprop_filter_ranges.mat'];
end
if(nargout < 2)
    imgBaseDir = pwd;
end
if(nargin < 3)
    manual = false;
end

%% Check that files exist
assert(manual || exist(rangesFile,'file') == 2,[rangesFile ' does not exist.']);
assert(exist(imgBaseDir,'dir') == 7,[rangesFile ' does not exist.']);

%% Load ranges data
% This file hould contain the structs native_regionprop_ranges and 
% synex_regionprop_ranges. The ranges would normally be chosen emperically
% by the user via a GUI. This contains sample values for analyzing the
% sample images.
if(~manual)
    ranges = load(rangesFile);
end

%% Load images
for i=1:4
    native.phase{i} = imread([imgBaseDir filesep 'native' filesep 'native_phase_' num2str(i) '.tif']);
    native.cfp{i} = imread([imgBaseDir filesep 'native' filesep 'native_cfp_' num2str(i) '.tif']);
end

for i=1:5
    synex.phase{i} = imread([imgBaseDir filesep 'synex' filesep 'synex_phase_' num2str(i) '.tif']);
    synex.cfp{i} = imread([imgBaseDir filesep 'synex' filesep 'synex_cfp_' num2str(i) '.tif']);
end

%% Do edge detection and retrieve connected components
native_cc = cellfun(@edgeDetect,native.phase,'UniformOutput',false);
synex_cc = cellfun(@edgeDetect,synex.phase,'UniformOutput',false);

%% Acquire region properties
region_properties = {'MeanIntensity','MinIntensity','MaxIntensity','Area','Solidity','Extent','EulerNumber','Perimeter','Eccentricity'};
native_rp = cellfun(@(I,CC) regionprops(CC,I,region_properties),native.phase,native_cc,'UniformOutput',false);
synex_rp = cellfun(@(I,CC) regionprops(CC,I,region_properties),synex.phase,synex_cc,'UniformOutput',false);

%% Filter connected component objects

% If we were manually doing this, then use this GUI interface
if(manual)
    native.f= cellfun(@(I,cc,rp) burstfiltering_s(I,'cc',cc','rp',rp),native.phase,native_cc,native_rp,'UniformOutput',false);
    synex.f= cellfun(@(I,cc,rp) burstfiltering_s(I,'cc',cc','rp',rp),synex.phase,native_cc,native_rp,'UniformOutput',false);
    
    native_regionprop_ranges = cellfun(@(f) f.ranges,native.f,'UniformOutput',false);
    native_regionprop_ranges = [native_regionprop_ranges{:}];
    synex_regionprop_ranges = cellfun(@(f) f.ranges,synex.f,'UniformOutput',false);
    synex_regionprop_ranges = [synex_regionprop_ranges{:}];
    
    save(rangesFile,'native_regionprop_ranges','synex_regionprop_ranges');
else
    [native_cc,native_rp] = cellfun(@applyFilter,native_cc,native_rp,num2cell(ranges.native_regionprop_ranges),'UniformOutput',false);
    [synex_cc,synex_rp] = cellfun(@applyFilter,synex_cc,synex_rp,num2cell(ranges.synex_regionprop_ranges),'UniformOutput',false);

    % Pack into structure, leaving out the output field which would have been
    % the handle to the burstfiltering figure
    native.f = cellfun(@(ranges,cc,rp) struct('ranges',ranges,'cc',cc,'rp',rp), ...
        num2cell(ranges.native_regionprop_ranges),native_cc,native_rp,'UniformOutput',false);

    synex.f = cellfun(@(ranges,cc,rp) struct('ranges',ranges,'cc',cc,'rp',rp), ...
        num2cell(ranges.synex_regionprop_ranges),synex_cc,synex_rp,'UniformOutput',false);
end

%% Generate masks from remaining connected components

native.mask = cellfun(@(f) labelmatrix(f.cc) > 0,native.f,'UniformOutput',false);
synex.mask = cellfun(@(f) labelmatrix(f.cc) > 0,synex.f,'UniformOutput',false);

%% Rescale masks since cfp channel is binned

native.mask = cellfun(@(mask) imresize(mask,0.5),native.mask,'UniformOutput',false);
synex.mask = cellfun(@(mask) imresize(mask,0.5),synex.mask,'UniformOutput',false);

%% Use mask to grab pixels from cfp channel

native.px = cellfun(@(cfp,mask) cfp(mask(:)),native.cfp,native.mask,'UniformOutput',false);
synex.px = cellfun(@(cfp,mask) cfp(mask(:)),synex.cfp,synex.mask,'UniformOutput',false);

end

