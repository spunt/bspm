function k = bspm_cluster_fwe(L2,L1,alpha,u,mask)
% BSPM_CLUSTER_FWE Calculate Extent for Cluster FWE Correction 
%
% USAGE: k = bspm_cluster_fwe(L1,L2,alpha,u,mask)
%
% INPUTS
%   L2: path to statistic image for level 2 model
%   L1: path to SPM.mat for level 1 models
%   alpha: corrected alpha level (default = .05)
%   u: uncorrected cluster-forming threshold (default = .001)
%   mask: search volume (default = all non-zero voxels in L2 statistic image)
%
% OUTPUT
%   k: corrected cluster size (spatial extent) threshold
%

% ---------------------------- Copyright (C) 2014 ----------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 5, mask = []; end
if nargin < 4, u = .001; end
if nargin < 3 || isempty(alpha), alpha = .05; end
if nargin < 2, L1 = []; end
if nargin < 1, disp('USAGE: k = bspm_cluster_fwe(level1mat,level2im,alpha,u,mask)'); return; end

%% DF %%
if iscell(L2), L2 = char(L2); end
img.hdr = spm_vol(L2);
tmp = img.hdr.descrip; i1 = find(tmp=='['); i2 = find(tmp==']');
df = str2num(tmp(i1(1)+1:i2(1)-1));
df = [1 df];
M = img.hdr.mat;

%% LEVEL 1 %%
if isempty(L1)
    [path name ext] = fileparts(L2);
    if isempty(path), path = pwd; end
    if exist([path filesep 'SPM.mat'],'file') 
        load([path filesep 'SPM.mat']);
        L1con = SPM.xY.P; 
        for c = 1:length(L1con)
            [p n e] = fileparts(L1con{c});
            load([p filesep 'SPM.mat']);
            subFWHM(c,:) = SPM.xVol.FWHM;
        end
    elseif exist([path filesep 'I.mat'],'file')
        load([path filesep 'I.mat']);
        L1con = I.F.IN.Scans; 
        for c = 1:length(L1con)
            [p n e] = fileparts(L1con{c});
            load([p filesep 'SPM.mat']);
            subFWHM(c,:) = SPM.xVol.FWHM;
        end
    else disp('Cannot find SPM.mat file! Quitting...'); return; end
else
    if ischar(L1), L1 = cellstr(L1); end
    for i = 1:length(L1);
        load(L1{i});
        subFWHM(i,:) = SPM.xVol.FWHM;
    end
end
FWHM = nanmean(subFWHM);
STAT = 'T';
range = 2:500;


%% MASK %%
if isempty(mask), 
    [path name ext] = fileparts(L2);
    if isempty(path), path = pwd; end
    if exist([path filesep 'mask.img'],'file'), mask.hdr = spm_vol([path filesep 'mask.img']); 
    elseif exist([path filesep 'mask.nii'],'file'), mask.hdr = spm_vol([path filesep 'mask.nii']);
    else disp('Cannot find mask file! Quitting...'); return; end
elseif iscell(mask),
    mask.hdr = spm_vol(char(mask)); 
else
    mask.hdr = spm_vol(mask);
end
mask.data = spm_read_vols(mask.hdr);
S = sum(mask.data(:)>0);
R = spm_resels_vol(mask.hdr,FWHM)';
v2r  = 1/prod(FWHM(~isinf(FWHM)));

%% CALCULATE %%
if u <= 1, u = spm_u(u,df,STAT); end
Pc = 1;
for k = range
    Pc = spm_P(1,k*v2r,u,df,STAT,R,1,S);
    if Pc <= alpha, 
      break 
    end
end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
