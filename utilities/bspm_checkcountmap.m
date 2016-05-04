function bspm_checkcountmap(countmap)
% BSPM_CHECKCOUNTMAP
% 
% USAGE: bspm_checkcountmap(countmap)
%
% ARGUMENTS
%   countmap = countmap of valid voxels 
%

% --------- Copyright (C) 2014 ---------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<1, mfile_showhelp; return; end
if iscell(countmap), countmap = char(countmap); end

% masks
maskdir = fullfile(getenv('HOME'), 'Github', 'bspm', 'imagedata', 'masks');
maskidx{1} = get_mask_idx(countmap);
maskidx{2} = get_mask_idx([maskdir filesep 'mask_amygdala_LR_edit.nii']);
maskidx{3} = get_mask_idx([maskdir filesep 'mask_vPFC.nii']);
maskidx{4} = get_mask_idx([maskdir filesep 'mask_vTP.nii']);
masknames = {'Full Mask' 'Amygdala' 'Ventral PFC' 'Ventral TP'};

% load count map
h = spm_vol(countmap);
cmap = spm_read_vols(h);
cmap = cmap(:);

% plot
n = max(cmap);
nplot = 1:n;
figure('Color','white');
count = 0;
for m = 1:length(maskidx)
    cm = cmap(maskidx{m});
    for i = nplot
        nvalid(i) = length(find(cm>=nplot(i)));
    end
    count = count + 1;
    axisarg = [0 nplot(end) min(nvalid)-(max(nvalid)/100) max(nvalid)+(max(nvalid)/100)];
    subplot(2,2,count);
    plot(nvalid,'LineWidth',1);
    set(get(gca,'YLabel'), 'String', 'N Valid Voxels', 'FontSize', 13);
    set(get(gca,'XLabel'), 'String', 'N Subjects', 'FontSize', 13);
    set(get(gca,'Title'), 'String', masknames{count}, 'FontSize', 15);
%     set(get(gca,'Axis'), 
    axis(axisarg);
end
function idx = get_mask_idx(maskfile)
hdr = spm_vol(maskfile);
img = spm_read_vols(hdr);
img = img(:);
idx = find(img);
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
