function out = bspm_harvest_betaseries(analysisdirs, patterns, rois)
% BSPM_HARVEST_BETA
%
%   USAGE: [data betanam roinam] = bspm_harvest_contrast(analysisdirs, rois, conidx, plotflag)
%
%   ARGUMENTS
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

% -------------------------------------- Copyright (C) 2014 --------------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<1, error('USAGE: bspm_harvest_beta(analysisdirs, rois, conidx, plotflag)'); end
if ischar(analysisdirs), analysisdirs = cellstr(analysisdirs); end
if ischar(rois), rois = cellstr(rois); end
nroi = length(rois);
nsub = length(analysisdirs);
roidata = bspm_read_vol(rois, 'reshape');
[~,roiname,~] = cellfun(@fileparts, rois, 'Unif', false);
subname = cellstrunique(analysisdirs);
for i = 1:nsub
    
    load(fullfile(analysisdirs{i}, 'SPM.mat'));
    allbname = {SPM.Vbeta.descrip}';
    allfname = {SPM.Vbeta.fname}';
    bidx = cellismember(allbname, patterns);
    % item names
    bname = regexp(allbname(bidx), 'Sn\(.\).+', 'match');
    bname = clean_text([bname{:,1}]', {'Sn\(\d\)' '\*bf\(\d\)'});
    % data
    fname = strcat(analysisdirs{i},filesep,allfname(bidx));
    cdata = bspm_read_vol(fname, 'reshape');
    subdata = zeros(length(fname), nroi);
    for r = 1:nroi
        subdata(:,r) = nanmean(cdata(find(roidata(:,r)),:))';
    end
    out{i}.subname = subname{i};
    out{i}.betaname = bname;
    out{i}.roiname = roiname; 
    out{i}.data = subdata;
    
end
    
    
    
   
 
 
 
 
