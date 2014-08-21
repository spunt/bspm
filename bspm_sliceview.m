function o = bspm_sliceview(under,over,view,slices,thresh,labels,cbar)
% BSPM_SLICEVIEW
%
% USAGE: o = bspm_sliceview(under,over,view,slices,scale)
%
%   ARGUMENTS
%       under: underlay
%       over: overlay
%       view: 'axial', 'coronal', or 'sagittal'
%       slices: vector of slices to display
%       scale: scale for overlay
%

% -------------------------- Copyright (C) 2014 --------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<7, cbar = 1; end
if nargin<6, labels = 0; end
if nargin<5, thresh = [.001 20]; end
if ischar(under), under = cellstr(under); end
if ischar(over), over = cellstr(over); end

% setup figure
% close all
set(0,'units','pixels');
pos = get(0, 'screensize');
pos(1:2) = 100;
pos(3:4) = pos(3:4) - 100;

% tmp filename 
tmpfile = [pwd filesep 'tmpfile.nii'];

for i = 1:length(over)

    % get overlay information
    data = bspm_threshold_image(over{i},thresh(1),thresh(2),0,tmpfile);
    [data hdr info] = bspm_read_vol(tmpfile);

    % initialize figure nad make slover object
    figure('Color','black','Position',pos);
    o = slover(char([under; cellstr(tmpfile)]));
    o.figure = [];
    o.slices = slices; 
    o.xslices = [];
    o.transform = lower(view);
    if ~labels, o.labels = 'none'; end
    if cbar, o.cbar = 2; else o.cbar = []; end
    o.refreshf = 0;
    o.resurrectf = 1;
    
    o.img(1).hold = 0;
    o.img(1).prop = 1;
    o.img(2).hold = 0;
    o.img(2).prop = 1;
    o.img(2).type = 'split';
    o.img(2).range = [info.min info.max]';
    o.img(2).cmap = hot(size(o.img(1).cmap,1));
    o.img(2).outofrange{1} = 0;
    
    o = paint(o);

end




 
 
 
 
