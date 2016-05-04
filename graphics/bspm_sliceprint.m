function o = bspm_sliceprint(under,over,view,slices,thresh,basename,labels,cbar,cmap)
% BSPM_SLICEPRINT
%
% USAGE: o = bspm_sliceview(under,over,view,slices,[thresh],[basename],[labels],[cbar],[cmap])
%
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
if nargin<9, cmap = []; end
if nargin<8, cbar = 0; end
if nargin<7, labels = 0; end
if nargin<6, basename = []; end
if nargin<5, thresh = [.001 20]; end
if ischar(under), under = cellstr(under); end
if ischar(over), over = cellstr(over); end

%% base name
if isempty(basename)
    [p n e] = fileparts(over{1});
    basename = upper(n);
end

%% figure
set(0,'units','pixels');
pos = get(0, 'screensize');
pos(1:2) = 100;
pos(3:4) = floor(pos(3:4)*.5);

%% tmp filename 
tmpfile = [pwd filesep 'tmpfile.nii'];
data = bspm_threshold_image(over{1},thresh(1),thresh(2),0,tmpfile);
tmp = unique(data(data>0));
if length(tmp)>10, cmapsize=64; else, cmapsize = length(tmp); end
[data hdr info] = bspm_read_vol(tmpfile);
%% colormap
if ischar(cmap)
    if length(tmp)>12, cmapsize=64; else cmapsize = length(tmp); end
end

[data hdr info] = bspm_read_vol(tmpfile);

for i = 1:length(slices)
    
    %% initialize figure and make slover object
    figure('color','white','position',pos);
    o = slover(char([under; cellstr(tmpfile)]));
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
    if ~labels, o.labels = 'none'; end
    if cbar, o.cbar = 2; else o.cbar = []; end
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
    o.img(2).range = [info.min info.max]';
    if isempty(cmap)
        o.img(2).cmap = hot(cmapsize);
    elseif ischar(cmap)
        o.img(2).cmap = brewermap(cmapsize,cmap);
    else
        o.img(2).cmap = cmap;
    end
    o.img(2).outofrange{1} = 0;
    currento = paint(o);
    tobj = findobj(gcf, 'type', 'text');
    set(tobj, 'units', 'data', 'fontsize', .15, 'string', regexprep(get(tobj, 'string'), ' ', '')); 
   
    %% save
    name = sprintf('%s_%02d_%s_%02d.jpg',basename,i,view,slices(i));
    bspm_save_figure(name, 'crop', 1); 
    close all
    
end

% %% colorbars
% set(0,'units','pixels');
% pos = get(0, 'screensize');
% pos(1:2) = 100;
% pos(3:4) = floor(pos(3:4)*.25);
% pos(4) = floor(pos(3)*.125);
% fs = floor(pos(4)*.5);
% figure('color','white','position',pos);
% imagesc(1:size(o.img(2).cmap,1)); 
% colormap(o.img(2).cmap);
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% set(gca,'xticklabel',[]);
% set(gca,'yticklabel',[]);
% % x = get(gca,'XLim');
% % y = get(gca,'YLim');
% % text(x(1),sum(y)/2,sprintf('%2.1f',info.min),'FontUnits','pixels','FontName','Arial', ...
% %     'HorizontalAlignment','Left','VerticalAlignment','Middle','FontSize',fs,'Color',[1 1 1]);
% % text(x(end),sum(y)/2,sprintf('%2.1f',info.max),'FontUnits','pixels','FontName','Arial', ...
% %     'HorizontalAlignment','Right','VerticalAlignment','Middle','FontSize',fs,'Color',[0 0 0]);
% outname = sprintf('%s_%2.2fto%2.2f_colorbar.png',basename,info.min,info.max);
% bspm_save_figure(outname);
% close all
% 



 
 
 
 
 
 
 
 
