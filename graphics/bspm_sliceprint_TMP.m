
under = {'/Users/bobspunt/Drive/Writing/Empirical/AVT/data/brain/meanT1_N=54.nii'};
over = {'/Users/bobspunt/Drive/Writing/Empirical/AVT/data/brain/nuisance_regressions/rt/spmT_0001.img'};
over2 = {'/Users/bobspunt/Drive/Writing/Empirical/AVT/data/brain/nuisance_regressions/rt/spmT_0002.img'};
slices = -6;
thresh = [.005 100];
view = 'sagittal';
labels = 1;
cbar = 0;

%% colormaps
blue = brewermap(100,'Blues');
red = brewermap(100,'Reds');


%% figure
set(0,'units','pixels');
pos = get(0, 'screensize');
pos(1:2) = 100;
pos(3:4) = pos(3:4) - 100;

%% tmp filename 
posfile = [pwd filesep 'posfile.nii'];
data = bspm_threshold_image(over{1},thresh(1),thresh(2),0,posfile);
[data hdr info] = bspm_read_vol(posfile);
negfile = [pwd filesep 'negfile.nii'];
data = bspm_threshold_image(over2{1},thresh(1),thresh(2),0,negfile);
[dataneg hdrneg infoneg] = bspm_read_vol(negfile);

for i = 1:length(slices)
    
    %% initialize figure and make slover object
    figure('color','white','position',pos);
    o = slover(char([under; cellstr(posfile); cellstr(negfile)]));
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
    
%     close all
%     clear hf name im imleft

end




 
 
 
 
 
 
 
 
 
 
 
 
