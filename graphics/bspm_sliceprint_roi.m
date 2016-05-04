function o = bspm_sliceprint_roi(under,over,view,slices,basename)
% BSPM_SLICEPRINT
%
% USAGE: o = bspm_sliceprint_roi(under,over,view,slices,basename)
%   ARGUMENTS
%       under: underlay
%       over: overlay
%       view: 'axial', 'coronal', or 'sagittal'
%       slices: vector of slices to display
%       thresh: [u k], e.g., [.001 20]
%       basename: for output filenames
%       labels: include labels (default = 0)
%       cbar: include color bar (default = 1)
%       cmap: colormap (default = hot)
%

% -------------------------------------- Copyright (C) 2014 --------------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<4, mfile_showhelp; return; end
if nargin<5, basename = []; end
if ischar(under), under = cellstr(under); end
if ischar(over), over = cellstr(over); end
if isempty(basename)
    [p, n, e] = fileparts(over{1});
    basename = upper(n);
end

%% figure
set(0,'units','pixels');
pos = get(0, 'screensize');
pos(1:2) = 100;
pos(3:4) = floor(pos(3:4)*.5);
    
%% initialize figure and make slover object
figure('color','white','position',pos);
o = slover(char([under; over]));
o.figure = [];
o.slices = slices; 
o.xslices = [];
view = lower(view);
if strcmp(view,'coronal');
    o.transform = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
    o.slicedef = [-78 1 78; -50 1 85];
elseif strcmp(view,'sagittal');
    o.transform = [0 -1 0 0; 0 0 1 0; 1 0 0 0; 0 0 0 1];
    o.slicedef = [-76 1 112; -50 1 85];
end
o.cbar = []; 
o.refreshf = 0;
o.resurrectf = 1;
o.area.units = 'normalized';
o.area.position = [0 0 1 1];
o.area.halign = 'left';
o.area.valign = 'top';
o.img(1).hold = 0;
o.img(1).prop = 1;
o.img(2).hold = 0;
o.img(2).prop = 1;
o.img(2).type = 'split';
o.img(2).range = [.5 1]; 
o.img(2).outofrange{1} = 0;
currento = paint(o);
tobj = findobj(gcf, 'type', 'text');
set(tobj, 'units', 'data', 'fontsize', .15, 'string', regexprep(get(tobj, 'string'), ' ', '')); 
name = sprintf('%s_%s_%02d.jpg',basename,view,slices);
bspm_save_figure(name, 'crop', 1); 
end
 
 
 
 
