function bspm_sliceprint(under,over,view,slices)
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

% ------------------- Copyright (C) 2014 -------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if ischar(under), under = cellstr(under); end
if ischar(over), over = cellstr(over); end
% if size(scale,1)==1, scale = scale'; end
set(0,'units','pixels');
pos = get(0, 'screensize');
w = pos(3) - pos(1); h = pos(4) - pos(2);
x = round(w*.5);
y = round(x/(4/3));
xgap = floor((w-x));
ygap = floor((h-y));

% determine scale automatically
[d hdr i] = bspm_read_vol(over);
scale = [i.min i.max]';

for i = 1:length(slices)
    
    figure('color', 'white','position',[xgap ygap (w-xgap) (h-ygap)]);
    o = slover(char([under; over]));
    o.figure = [];
    o.slices = slices(i);
    o.xslices = [];
    view = lower(view);
    if strcmp(view,'coronal');
        o.transform = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
        o.slicedef = [-78 1 78; -50 1 85];
    elseif strcmp(view,'sagittal');
        o.transform = [0 -1 0 0; 0 0 1 0; 1 0 0 0; 0 0 0 1];
        o.slicedef = [-76 1 112; -50 1 85];
    end
    o.refreshf = 0;
    o.resurrectf = 1;
    o.img(1).hold = 0;
    o.img(1).prop = 1;
    o.img(2).hold = 0;
    o.img(2).prop = 1;
    o.img(2).type = 'split';
    o.img(2).range = scale;
    o.img(2).cmap = hot(size(o.img(1).cmap,1));
    o.img(2).outofrange{1} = 0;
    o.cbar = []; % no color bar if empty (2 = colorbar)
%     o.xslices = 2;
    o.labels = 'none'; 
    paint(o);
    
    % get and correct
    hf = gcf;
    im = getframe(hf);
    im = im.cdata;
    dim = size(im)/2;
    figidx = find(im(:,end,1)<100);
    imleft = im(figidx,1:dim(2)-1,:);

    % save
    name = sprintf('%03d_%s_%d.jpg',i,view,slices(i));
    imwrite(imleft, name, 'jpg');
    close all
    clear hf name im imleft
    
    
end




 
 
 
 
