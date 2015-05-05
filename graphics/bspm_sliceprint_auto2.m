function o = bspm_sliceprint_auto2(under,over,view,thresh,basename,labels,cbar)
% BSPM_SLICEPRINT
%
% USAGE: o = bspm_sliceview(under,over,view,slices,[thresh],[basename],[labels],[cbar])
%
%   ARGUMENTS
%       under: underlay
%       over: overlay
%       view: 'axial', 'coronal', or 'sagittal'
%       slices: vector of slices to display
%       thresh: [u k sep], e.g., [.001 20 18]
%       basename: for output filenames
%       labels: include labels (default = 0)
%       cbar: include color bar (default = 1)
%

% ---------------------------------- Copyright (C) 2014 ----------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<7, cbar = 0; end
if nargin<6, labels = 1; end
if nargin<7, basename = []; end
if nargin<4, thresh = [.001 20 18]; end
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
posfile = [pwd filesep 'posfile.nii'];
data = bspm_threshold_image(over{1},thresh(1),thresh(2),0,posfile);
[data hdr info] = bspm_read_vol(posfile);

%% negative
[tmpdata tmphdr tmpinfo] = bspm_read_vol(over{1});
negfile = [pwd filesep 'negfile.nii'];
negfile2 = [pwd filesep 'negfile_thresh.nii'];
tmphdr.fname = negfile;
tmpdata = tmpdata*-1;
spm_write_vol(tmphdr,tmpdata);
data = bspm_threshold_image(negfile,thresh(1),thresh(2),0,negfile2);
[dataneg hdrneg infoneg] = bspm_read_vol(negfile2);

%% get slices
pospeaks = bspm_get_peaks(over{1},thresh(1),thresh(2),thresh(3));
negpeaks = bspm_get_peaks(negfile,thresh(1),thresh(2),thresh(3));
peaks = [pospeaks.coords; negpeaks.coords];

if strcmp(view,'coronal');
    slices = peaks(:,2);
elseif strcmp(view,'sagittal');
    slices = peaks(:,1);
else
    slices = peaks(:,3);
end
slices = unique(slices);

for i = 1:length(slices)
    
    %% initialize figure and make slover object
    figure('color','white','position',pos);
    o = slover(char([under; cellstr(posfile); cellstr(negfile2)]));
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
    o.img(2).cmap = brewermap(size(o.img(1).cmap,1),'Reds');
    o.img(2).outofrange{1} = 0;
    o.img(3).hold = 0;
    o.img(3).prop = 1;
    o.img(3).type = 'split';
    o.img(3).range = [infoneg.min infoneg.max]';
    o.img(3).cmap = brewermap(size(o.img(1).cmap,1),'Blues');
    o.img(3).outofrange{1} = 0;
    paint(o);
    
    %% get and correct
    hf = gcf;
    im = getframe(hf);
    im = im.cdata;
    im2 = rgb2gray(im);
    hidx = find(nanmean(im2)>200);
    vidx = find(nanmean(im2,2)>200); 
    im(vidx,:,:) = [];
    im(:,hidx,:) = [];


    %% save
    name = sprintf('%s_%02d_%s_%02d.jpg',basename,i,view,slices(i));
    imwrite(im, name, 'jpg');
    close all
    jheapcl

end
%% colorbars
set(0,'units','pixels');
pos = get(0, 'screensize');
pos(1:2) = 100;
pos(3:4) = floor(pos(3:4)*.25);
pos(4) = floor(pos(3)*.125);
fs = floor(pos(4)*.5);

figure('color','white','position',pos);
h1 = imagesc(1:size(o.img(2).cmap,1)); 
colormap(o.img(2).cmap);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
x = get(gca,'XLim');
y = get(gca,'YLim');
text(x(1),sum(y)/2,sprintf('%2.1f',info.min),'FontUnits','pixels','FontName','Arial', ...
    'HorizontalAlignment','Left','VerticalAlignment','Middle','FontSize',fs,'Color',[0 0 0]);
text(x(end),sum(y)/2,sprintf('%2.1f',info.max),'FontUnits','pixels','FontName','Arial', ...
    'HorizontalAlignment','Right','VerticalAlignment','Middle','FontSize',fs,'Color',[1 1 1]);
bob_save_figure('colorbar1.png');
close all

figure('color','white','position',pos);
h2 = imagesc(1:size(o.img(3).cmap,1)); 
colormap(o.img(3).cmap);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
x = get(gca,'XLim');
y = get(gca,'YLim');
text(x(1),sum(y)/2,sprintf('%2.1f',infoneg.min),'FontUnits','pixels','FontName','Arial', ...
    'HorizontalAlignment','Left','VerticalAlignment','Middle','FontSize',fs,'Color',[0 0 0]);
text(x(end),sum(y)/2,sprintf('%2.1f',infoneg.max),'FontUnits','pixels','FontName','Arial', ...
    'HorizontalAlignment','Right','VerticalAlignment','Middle','FontSize',fs,'Color',[1 1 1]);
bob_save_figure('colorbar2.png');
close all





 
 
 
 
