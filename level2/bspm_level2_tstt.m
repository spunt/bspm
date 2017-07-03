function matlabbatch = bspm_level2_tstt(cons1, cons2, varargin)
% BSPM_LEVEL2_TSTT
%
%   ARGUMENTS:
%       cons1: contrast images from level 1 for group 1
%       cons2: contrast images from level 1 for group 2
%       grouplabels: labels for group 1 and group 2 (e.g., {'group1' 'group2'})
%       covariates: specified as follows:
%           covariates(n).name = 'Covariate Name';
%           covariates(n).values = [1...n];
%

% -------------------------- Copyright (C) 2014 --------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
def =   { ...
        'outdir',       [],                     ...
        'grouplabels',  {'Group1' 'Group2'},    ...
        'tag',          [],                     ...
        'implicitmask',  1,                     ...
        'omitpat',      {''},   ...
        'mask',         '',                     ...
        'pctgroup',     [],                     ...
        'viewit',        0,                     ...
        'negativecon',   0,                     ...
        'nan2zero',      0,                     ...
        'covariates',   []                      ...
        };
vals = setargs(def, varargin);
if nargin<2, mfile_showhelp; fprintf('\t= DEFAULT SETTINGS =\n'); disp(vals); return; end
fprintf('\n\t= CURRENT SETTINGS =\n'); disp(vals);
if ischar(cons1), cons1 = cellstr(cons1); end
if size(cons1,1)==1, cons1 = cons1'; end
if ischar(cons2), cons2 = cellstr(cons2); end
if size(cons2,1)==1, cons2 = cons2'; end
if iscell(outdir), outdir = char(outdir); end
if iscell(mask), mask = char(mask); end

rmidx1 = cellismember(cons1, omitpat); cons1(rmidx1) = [];
rmidx2 = cellismember(cons2, omitpat); cons2(rmidx2) = [];
if ~isempty(mask)
    [~,mtag] = fileparts(mask);
else
    mtag = 'noexplicitmask';
end

% | OUTPUT DIRECTORY
if isempty(outdir)

    % | Contrast Name
    cname = char(bspm_con2name(cons1{1}));

    % | Analysis Name
    [p, level1name] = fileparts(fileparts(cons1{1}));
    gadir           = fullfile(parentpath(cons1), '_groupstats_', level1name);
    if tag
        gasubdir    = fullfile(gadir, sprintf('TSTT_%s_%s_N%d_vs_%s_N%d_%s_%s', tag, grouplabels{1}, length(cons1), grouplabels{2}, length(cons2), mtag, bspm_timestamp(1)));
    else
        gasubdir    = fullfile(gadir, sprintf('TSTT_%s_N%d_vs_%s_N%d_%s_%s', grouplabels{1}, length(cons1), grouplabels{2}, length(cons2), mtag, bspm_timestamp(1)));
    end
    outdir          = fullfile(gasubdir, cname);

    % | Make Directories
    if ~isdir(gadir), mkdir(gadir); end
    if ~isdir(gasubdir), mkdir(gasubdir); end

end
if ~isdir(outdir), mkdir(outdir); end

% | NAN2ZERO (IF APPLICABLE)
if nan2zero
    bspm_batch_imcalc([cons1; cons2], '', 'nan2zero');
end

% | PCTGROUP (IF APPLICABLE)
if ~isempty(pctgroup)
    cons = [cons1; cons2];
    if ~isempty(mask)
        [d,h] = bspm_read_vol(cons, 'mask', mask);
    else
        [d,h] = bspm_read_vol(cons);
    end
    d(isnan(d))               = 0;
    m                         = double(sum(d~=0, 4))/length(cons);
    if pctgroup > 1, pctgroup = pctgroup/100; end
    m(m<pctgroup)             = 0;
    mask                      = fullfile(outdir, sprintf('Mask_PctGroup%d_%s.nii', round(pctgroup*100), mtag));
    hdr                       = h(1);
    hdr.fname                 = mask;
    hdr.descrip               = sprintf('Valid Voxels Mask - PercentGroup=%d  - %s', round(pctgroup*100), mtag);
    spm_write_vol(hdr, m);
end

% | FIX END OF IMAGE FILENAMES
cons1 = strcat(cons1, ',1');
cons2 = strcat(cons2, ',1');
if ~isempty(mask), mask = strcat(mask, ',1'); end

% | FACTORIAL DESIGN SPECIFICATION
matlabbatch{1}.spm.stats.factorial_design.dir{1}                 = outdir;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1          = cellstr(cons1);
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2          = cellstr(cons2);
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept            = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance        = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca           = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova          = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none     = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im             = implicitmask;
matlabbatch{1}.spm.stats.factorial_design.masking.em{1}          = mask;
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit         = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm        = 1;

% | COVARIATES
if ~isempty(covariates)
    ncov = length(covariates);
    for i = 1:ncov
        matlabbatch{1}.spm.stats.factorial_design.cov(i).c      = covariates(i).values;
        matlabbatch{1}.spm.stats.factorial_design.cov(i).cname  = covariates(i).name;
        matlabbatch{1}.spm.stats.factorial_design.cov(i).iCFI   = 1;
        matlabbatch{1}.spm.stats.factorial_design.cov(i).iCC    = 1;
    end
end

% | ESTIMATE
matlabbatch{2}.spm.stats.fmri_est.spmmat{1}        = fullfile(outdir,'SPM.mat');
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% | CONTRASTS
matlabbatch{3}.spm.stats.con.spmmat{1}               = fullfile(outdir,'SPM.mat');
matlabbatch{3}.spm.stats.con.delete                  = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name    = [grouplabels{1} '_-_' grouplabels{2}];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec  = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
if negativecon
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name    = [grouplabels{2} '_-_' grouplabels{1}];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec  = [-1 1];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
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
