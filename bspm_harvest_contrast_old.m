function [data connam roinam] = bspm_harvest_contrast(analysisdirs, rois, conidx, plotflag)
% BSPM_HARVEST_CONTRAST
%
%   USAGE: [data connam roinam] = bspm_harvest_contrast(analysisdirs, rois, conidx, plotflag)
%
%   ARGUMENTS
%       analysisdirs: analysis directory containing contrast images
%       rois: paths to region of interest images from which to extract data
%       conidx: option to harvest a subset of contrasts indexed here
%       plotflag: option to plot bargraph of contrasts for each ROI (default = 0)
%
%   OUTPUTS
%       data: a cell array of the contrast values (each cell contains data
%       matrix for a different roi corresponding to order of roinam;
%       subjects are rows and contrasts are columns)
%       connam: a cell array of contrast names
%       roinam: a cell array of ROI names
%       
% Created March 29, 2013 - Bob Spunt

% ------------------------------------- Copyright (C) 2014 -------------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<4, plotflag = 0; end
if nargin<3, conidx = []; end
if nargin<2, error('USAGE: bspm_harvest_contrast(analysisdirs, rois, conidx, plotflag)'); end

% make sure image names are character arrays
% ------------------------------------------------------
if ischar(analysisdirs), analysisdirs = cellstr(analysisdirs); end
if ischar(rois), rois = cellstr(rois); end

% get indices for ROIs
% ------------------------------------------------------
nroi = length(rois);
for i = 1:nroi
    
    roi = rois{i};
    [path roinam{i} e] = fileparts(roi);
    hdr = spm_vol(roi); img = spm_read_vols(hdr);
    roiIDX{i} = find(img);
    
end

% loop over analyses
% ------------------------------------------------------
nanalysis = length(analysisdirs);

for i = 1:nanalysis

    spmmat = [analysisdirs{i} filesep 'SPM.mat'];
    tmp = load(spmmat);
    connam = {tmp.SPM.xCon.name};
    stat = {tmp.SPM.xCon.STAT};
    if ~isempty(conidx)
        connam = connam(conidx);
        stat = stat(conidx);
    end
    docon = cellstrfind(stat,'T');
    connam = connam(docon);
    connam = regexprep(connam,' - All Sessions', '');

    for r = 1:nroi
        
        for c = docon
            
            conimg = [analysisdirs{i} filesep 'con_' sprintf('%04d',c) '.img'];
            hdr = spm_vol(conimg); img = spm_read_vols(hdr);
            data{r}(i,c) = nanmean(img(roiIDX{r}));
                   
        end
        
    end
    
end

if plotflag
    
    plotdata = [];
    groupidx = [];
    for d = 1:length(data)
        
        pos = size(plotdata,2);
        plotdata = [plotdata data{d}];
        groupidx(d,:) = pos+1:pos+size(data{d},2);
        
    end
    fontsize = [10 12];
    label_y = '% Signal Change';
    label_x = 'Region of Interest';
    label_groups = roinam;
    label_title = '';
    label_legend = regexprep(connam,'_',' ');
    
    % plot!
    figure('Color', 'white')
    bob_bargraph(plotdata, groupidx, label_title, label_x, label_y, label_groups, label_legend, fontsize)
    
end
        










    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
 
 
 
 
