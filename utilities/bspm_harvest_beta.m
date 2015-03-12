function [data, condname, roiname] = bspm_harvest_beta(analysisdirs, roifiles, betaidx)
% BSPM_HARVEST_BETA
%
%   USAGE: [data condname roiname] = bspm_harvest_beta(analysisdirs, rois, betaidx)
%
%   ARGUMENTSå
%       analysisdirs: analysis directory containing contrast images
%       rois: paths to region of interest images from which to extract data
%
%   OUTPUTS
%       data: a cell array of the contrast values (each cell contains data
%       matrix for a different roi corresponding to order of roinam;
%       subjects are rows and contrasts are columns)
%       beta: a cell array of beta names
%       roinam: a cell array of ROI names
%       
% Created March 29, 2013 - Bob Spunt

% ---------------------------------- Copyright (C) 2014 ----------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<2, error('USAGE: [data, condname, roiname] = bspm_harvest_beta(analysisdirs, roifiles)'); end
if nargin<3, betaidx = []; end

% make sure image names are character arrays
% ------------------------------------------------------
if ischar(analysisdirs), analysisdirs = cellstr(analysisdirs); end
if ischar(roifiles), roifiles = cellstr(roifiles); end

% get indices for ROIs
% ------------------------------------------------------
nroi = length(roifiles);
for i = 1:nroi
    
    roi = roifiles{i};
    [path roiname{i} e] = fileparts(roi);
    hdr = spm_vol(roi); img = spm_read_vols(hdr);
    roiIDX{i} = find(img);
    
end

% loop over analyses
% ------------------------------------------------------
nanalysis = length(analysisdirs);

for i = 1:nanalysis

    spmmat = [analysisdirs{i} filesep 'SPM.mat'];
    tmp = load(spmmat);
    condname = tmp.SPM.xX.name;
    betanam = {tmp.SPM.Vbeta.fname};
    % get only convolved regressors
    if isempty(betaidx)
        idx = cellstrfind(condname,'bf(1)');
        condname = condname(idx);
        betanam = betanam(idx);
    else
        idx = betaidx; 
    end
    % clean up names
    for r = 1:length(tmp.SPM.Sess)
        string = sprintf('Sn\\(%d\\) ', r);
        condname = regexprep(condname,string,'');
    end
    for r = 1:nroi
        
        for c = 1:length(betanam)
            
            conimg = [analysisdirs{i} filesep betanam{c}];
            hdr = spm_vol(conimg); img = spm_read_vols(hdr);
            data{r}(i,c) = nanmean(img(roiIDX{r}));
                   
        end
        
    end
    
end
end
 
 
 
 
