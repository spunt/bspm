function [data connam roinam] = bspm_harvest_contrast(analysisdirs, rois, conidx, conpat, plotflag)
% BSPM_HARVEST_CONTRAST
%
%   USAGE: [data connam roinam] = bspm_harvest_contrast(analysisdirs, rois, conidx, conpat, plotflag)
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

% ---------------------------------------- Copyright (C) 2014 ----------------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<5, plotflag = 0; end
if nargin<4, conpat = 'con*img'; end
if nargin<3, conidx = []; end
if nargin<2, error('USAGE: bspm_harvest_contrast(analysisdirs, rois, conidx, conpat, plotflag)'); end

if ischar(rois), rois = cellstr(rois); end
nroi = length(rois);
[p,roinam,e] = cellfun(@fileparts,rois,'UniformOutput',false);
roidata = bspm_read_vol_multi(rois,1);
nanalysis = length(analysisdirs);
for i = 1:nanalysis
    
    con = files([analysisdirs{i} filesep conpat]);
    if ~isempty(conidx), con = con(conidx); end
    ncon = length(con);
    h = spm_vol(char(con));
    descrip = {h.descrip};
    for c = 1:length(descrip)
        cn = descrip{c};
        idx = strfind(cn,': ');
        cn = cn(idx+2:end);
        cn = regexprep(cn, '- All Sessions','');
        connam{c} = strtrim(cn);
    end
    condata = bspm_read_vol_multi(con,1);
    for r = 1:nroi
        for c = 1:ncon
            data{r}(i,c) = nanmean(condata(roidata(:,r)==1,c)); 
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
    fontsizes = [10 12 12];
    labels.y = '% Signal Change';
    labels.x = 'Region of Interest';
    labels.groups = roinam;
    labels.title = '';
    labels.legend = regexprep(connam,'_',' ');
    labels.legend_title = 'Condition';
    
    % plot!
    figure('Color', 'white')
    bob_bargraph(plotdata, groupidx, labels, fontsizes)
% OPTIONAL ARGUMENTS
%   groupidx: each row indexes a different grouping
%   labels: a structure with the following fields (all optional)
%       title = plot title
%       x = x-axis label
%       y = y-axis label
%       groups = labels for groupings of bars
%       legend = labels for bars (in each group if applicable)
%       legend_title = label for the legend title
%   fontsizes: 1 x 3 vector with sizes for
%       Cell 1 = x and y axis labels
%       Cell 2 = legend and group labels
%       Cell 3 = plot title 
end
        










    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
 
 
 
 
