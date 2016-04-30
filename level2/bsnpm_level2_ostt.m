function matlabbatch = bsnpm_level2_ostt(cons, varargin)
% BSNPM_LEVEL2_OSTT
%
%   USAGE: matlabbatch = bsnpm_level2_ostt(cons, varargin)
%
%       cons: contrast images from level 1
%
%       covariates: specified as follows: 
%           covariates(n).name = 'Covariate Name';
%           covariates(n).values = [1...n]; 
%
%   OPTIONAL
%       outdir: allows custom output directory
%     implicit: yes (1) or no (0) for implicit masking
%         mask: specify explicit mask file if desired
%     pctgroup: specify some percent of all subjects. voxels present in at
%     least that percentage of subjects will be included. default is to not do this. 
%       viewit: will change to output directory and open BSPMVIEW at finish
%     nan2zero: will convert nans to zeros in the con images (important for
%     explicit masking)
%
% HELP FROM SNPM_CP.m
%-----------------------------------------------------------------------
%
%
% SingleSub: Two Sample T test; 2 conditions, replications
% SingleSub: Simple Regression (correlation); single covariate of interest
% MultiSub: One Sample T test on diffs/contrasts; 1 condition, 1 scan per subject
% MultiSub: Simple Regression (correlation); single covariate of interest, 1 scan per subject
% MultiSub: Paired T test; 2 conditions, 1 scan per condition
% MultiSub: Within Subject ANOVA; multiple scans/subject
% 2 Groups: Test diff of response; 2 conditions, 1 scan per condition
% 2 Groups: Two Sample T test; 1 scan per subject
% >2 Groups: Between Group ANOVA; 1 scan per subject
%
% Output File Descriptions:
%
%   XYZ.mat contains a 3 x N matrix of the x,y and z location of the
% voxels in SPMF in mm (usually referring the the standard anatomical
% space (Talairach and Tournoux 1988)} (0,0,0) corresponds to the
% centre of the voxel specified by ORIGIN in the *.hdr of the original
% and related data.
%
%   BETA.mat contains a p x S matrix of the p parameter estimates at
% each of the S voxels for the correct permutation.  These parameters 
% include all effects specified by the design matrix.
%
%   SnPMt.mat contains a 1 x S matrix of the statistic of interest (either
% t or pseudo-t if variance smoothing is used) supplied for all S voxels at
% locations XYZ.
%
%   SnPMucp.mat contains a 1 x S matrix of the nonparametric P values of
% the statistic of interest supplied for all S voxels at locations XYZ.
%
%   SnPM.mat contains a collection of strings and matrices that pertain 
% to the analysis.  In contrast to spm_spm's SPM.mat, most of the essential
% matrices are in the any of the matrices stored here in the CfgFile
% and hence are not duplicated here.   Included are the number of voxels
% analyzed (S) and the image and voxel dimensions [V].  See below
% for complete listing.
%
% snpm_cp writes out the following image files (for each image, there are
% two files: .img and .hdr files)  
%  
% beta_**** (from 0001 to p): p images of p parameter estimates at each
% voxel for the correct permutation. These p parameters include all
% effects specified by the design matrix. 
%
% ResMS: One image of residual mean square errors at each voxel. 
% 
% (SnPM, like SPM, only implements single tailed tests. In the following
% files, '+' or '-' correspond to 'positive' or 'negative' effects (as in
% snpm_pp.m). Here, '+' images are the images for large values,
% indicating evidence against the null hypothesis in favour of a positive
% alternative (activation, or positive slope in a covariate analysis))
% 
% snpmT+ & snpmT-: Images of the statistic of interest (either t or
% pseduo-t if variance smoothing is used), positive or negative. 
% The numbers (i.e. not NaN) saved in snpmT+ images are also saved in the
% SnPMt.mat file. 
% 
% lP+ & lP-: Images of -log10(uncorrected non-parametric P-values,
% positive or negative).
% 
% lP_FWE+ & lP_FWE-: Images of -log10(FWE-corrected non-parametric
% P-values, positive or negative). Here, FWE-corrected non-parametric
% P-values are the proportion of the permutation distribution for the
% maximal statistic which exceeds the statistic image at the voxel. 
%
% lP_FDR+ & lP_FDR-: Images of -log10(FDR-corrected non-parametric
% P-values, positive or negative). 
%
% The following is an example of matlab codes for reading in an image file. 
% P='.../.../beta_0001.img';
% V=spm_vol(P);
% Y=spm_read_vols(V);
% Y(~isfinite(Y))=[]; %delete NaN values from vector Y.
%
%
%-----------------------------------------------------------------------
%
% As an "engine", snpm_cp does not produce any graphics; if the SPM windows
% are open, a progress thermometer bar will be displayed.  
%
% If out-of-memory problems are encountered, the first line of defense is to
% run snpm_cp in a virgin matlab session with out first starting SPM.
%
%
% Variables saved in SnPM.mat
%=======================================================================
%
% S              Volume analyzed (in voxels)
% V              Volume handles (see spm_vol)
% df             Residual degrees of freedom of raw t-statistic
% MaxT           nPerm x 2 matrix of [max;min] t-statistics per perm
% ST_Ut          Threshold above which suprathreshold info was collected.
%                Voxel locations, t and perm are saved in SnPM_ST.mat for
%                t's greater than ST_Ut. ST_Ut=Inf if not saving STCdata
%
% s_SnPM_save    List of variables saved in SnPM.mat file
% CfgFile        SnPM config sile used (full pathname)
% s_SnPMcfg_save List of variables saved in SnPMcfg.mat file
% 
% Data structure of SnPM_ST.mat: suprathreshold stats (if collected)
%-----------------------------------------------------------------------
% 5xn matrix, each column containing:
%       [x, y, z, abs(T), perm]'
%       perm is negative if T was negative
%

% ---------------------------------- Copyright (C) 2014 ----------------------------------
%   Author: Bob Spunt
%   Affilitation: Caltech
%   Email: spunt@caltech.edu
%
%   $Revision Date: Aug_20_2014
def = { ...
        'clusterFormingThresh', .001,           ...
        'covariates',           [],             ...
        'implicitmask',         0,              ...
        'level2correct',        'cluster',      ... 
        'mask',                 '',             ...
        'nan2zero',             1               ...
        'nPerm',                10000,          ...
        'outdir',               [],             ...
        'pctgroup',             90,             ...
        'tag',                  [],             ...
        'varianceSmoothing',    [0 0 0],        ...
        'viewit',               0,              ...
        };
vals = setargs(def, varargin);
if nargin==0, mfile_showhelp; fprintf('\t= DEFAULT SETTINGS =\n'); disp(vals); return; end
fprintf('\n\t= CURRENT SETTINGS =\n'); disp(vals);
if length(varianceSmoothing)==1, varianceSmoothing = repmat(varianceSmoothing, 1, 3); end 
if ischar(cons), cons = cellstr(cons); end
if iscell(mask), mask = char(mask); end
if ~isempty(mask)
    [~,mtag] = fileparts(mask);
    mtag = upper(mtag);
else
    mtag = 'NOMASK'; 
end

% | NAN2ZERO (IF APPLICABLE)
if nan2zero, bspm_batch_imcalc(cons, '', 'nan2zero'); end

% | OUTPUT DIRECTORY
if isempty(outdir)
    
    % | Contrast Name
    cname = char(bspm_con2name(cons{1}));
    
    % | Analysis Name
    [p, level1name]  = fileparts(fileparts(cons{1})); 
    gadir       = fullfile(parentpath(cons), '_groupstats_', level1name);
    if isempty(pctgroup), pctgrouptag = 100; else pctgrouptag = pctgroup; end
    if tag
        gasubdir    = fullfile(gadir, sprintf('SnPM_%s_OSTT_%s_%dPERM_vFWHM%dx%dx%d_N%d_PCTIN%d_%s_%s', upper(level2correct), tag, nPerm, varianceSmoothing, length(cons), pctgrouptag, mtag, bspm_timestamp(1)));
    else
        gasubdir    = fullfile(gadir, sprintf('SnPM_%s_OSTT_%dPERM_vFWHM%dx%dx%d_N%d_PCTIN%d_%s_%s', upper(level2correct), nPerm, varianceSmoothing, length(cons), pctgrouptag, mtag, bspm_timestamp(1)));
    end
    outdir      = fullfile(gasubdir, cname);

    % | Make Directories
    if ~isdir(gadir), mkdir(gadir); end
    if ~isdir(gasubdir), mkdir(gasubdir); end
    
end
if ~isdir(outdir), mkdir(outdir); end

% | PCTGROUP (IF APPLICABLE)
if ~isempty(pctgroup)
    if ~isempty(mask)
        [d,h] = bspm_read_vol(cons, 'mask', mask);
    else
        [d,h] = bspm_read_vol(cons); 
    end
    d(isnan(d)) = 0; 
    m = double(sum(d~=0, 4))/length(cons);
    if pctgroup > 1, pctgroup = pctgroup/100; end
    m(m<pctgroup) = 0;
    mask = fullfile(outdir, sprintf('Mask_PctGroup%d_%s.nii', round(pctgroup*100), mtag));
    hdr = h(1); 
    hdr.fname = mask; 
    hdr.descrip = sprintf('Valid Voxels Mask - PercentGroup=%d  - %s', round(pctgroup*100), mtag); 
    spm_write_vol(hdr, m); 
end

% | FIX END OF IMAGE FILENAMES
cons = strcat(cons, ',1'); 
if ~isempty(mask), mask = strcat(mask, ',1'); end

% | DESIGN SPECIFICATION
matlabbatch{1}.spm.tools.snpm.des.OneSampT.DesignName = 'MultiSub: One Sample T test on diffs/contrasts';
matlabbatch{1}.spm.tools.snpm.des.OneSampT.DesignFile = 'snpm_bch_ui_OneSampT';
matlabbatch{1}.spm.tools.snpm.des.OneSampT.dir{1}     = outdir; 
matlabbatch{1}.spm.tools.snpm.des.OneSampT.P          = cellstr(cons);
matlabbatch{1}.spm.tools.snpm.des.OneSampT.nPerm      = nPerm;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.vFWHM      = varianceSmoothing;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.bVolm      = 1;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.ST.ST_U    = clusterFormingThresh;

% | MASKING & GLOBAL CALCULATIONS
matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.tm.tm_none            = 1;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.im                    = implicitmask;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.em{1}                 = mask;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalc.g_omit                = 1;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.gmsca.gmsca_no        = 1;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.glonorm               = 1;

% | COVARIATES
if ~isempty(covariates)
    ncov = length(covariates);
    for i = 1:ncov
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.cov(i).c      = covariates(i).values;
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.cov(i).cname  = covariates(i).name;
    end
end

% | ESTIMATION
matlabbatch{2}.spm.tools.snpm.cp.snpmcfg(1) = cellstr(fullfile(outdir,'SnPMcfg.mat'));

switch lower(level2correct)
    case 'cluster'
        % | INFERENCE (CLUSTER)
        matlabbatch{3}.spm.tools.snpm.inference.SnPMmat(1)                       = cellstr(fullfile(outdir,'SnPM.mat'));
        matlabbatch{3}.spm.tools.snpm.inference.Thr.Clus.ClusSize.CFth           = NaN;
        matlabbatch{3}.spm.tools.snpm.inference.Thr.Clus.ClusSize.ClusSig.FWEthC = 0.05;
        matlabbatch{3}.spm.tools.snpm.inference.Tsign                            = 1;
        matlabbatch{3}.spm.tools.snpm.inference.WriteFiltImg.name                = 'SnPM_ClusterLevel_Filtered';
        matlabbatch{3}.spm.tools.snpm.inference.Report                           = 'MIPtable'; 
    otherwise
        % | INFERENCE (VOXEL)
        matlabbatch{3}.spm.tools.snpm.inference.SnPMmat(1)                       = cellstr(fullfile(outdir,'SnPM.mat'));
        matlabbatch{3}.spm.tools.snpm.inference.Thr.Vox.VoxSig.FWEth             = 0.05;
        matlabbatch{3}.spm.tools.snpm.inference.Tsign                            = 1;
        matlabbatch{3}.spm.tools.snpm.inference.WriteFiltImg.name                = 'SnPM_VoxelLevel_Filtered';
        matlabbatch{3}.spm.tools.snpm.inference.Report                           = 'MIPtable';
end
% | RUN IF NO OUTPUT ARGS SPECIFIED
if nargout==0
    bspm_runbatch(matlabbatch);
    if viewit
        cd(outdir); 
        bspmview; 
    end
end

end



