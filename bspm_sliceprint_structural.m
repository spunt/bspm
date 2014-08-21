function o = bspm_sliceprint_structural(struct, view, slices, labels, name)
% BSPM_SLICEVIEW
%
% USAGE: o = bspm_sliceprint_structural(anat, view, slices, labels, name)
%
%   ARGUMENTS
%       struct: structural image
%       view: 'axial', 'coronal', or 'sagittal'
%       slices: vector of slices to display
%       labels: tag to include labels (default = 0)
%       name: optional name 
%

% ---------------------------- Copyright (C) 2014 ----------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<5, name = []; end
if nargin<4, labels = 0; end
if nargin<3, error('USAGE: o = bspm_sliceprint_structural(struct, view, slices, labels)'); end
if ischar(struct), anat = cellstr(struct); end

% set up figure
set(0,'units','pixels');
pos = get(0,'screensize');
pos(1:2) = 100;
pos(3:4) = pos(3:4) - 100;

for i = 1:length(struct)

    anat = struct{i};
    figure('Color','black','Position',pos);

    % make slover object
    o = slover(anat);
    o.figure = [];
    o.slices = slices; 
    o.xslices = [];
    o.transform = lower(view);
    if ~labels, o.labels = 'none'; end
    o.cbar = [];
    o.refreshf = 0;
    o.resurrectf = 1;

    % paint slover object and get figure
    o = paint(o);
    im = getframe(gcf);
    im = im.cdata;
    dim = size(im)/2;
    figidx = find(im(:,end,1)<100);
    im = im(figidx,:,:);

    % save
    if isempty(name)
        [p n e] = fileparts(anat);
        fn = sprintf('sliceview_%s_%s.png',view,n);
    else
        fn = name;
    end
    imwrite(im, fn, 'png');
    
    close all
    clear o im dim 
    
end
 
 
 
 
