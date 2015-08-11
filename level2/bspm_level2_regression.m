function matlabbatch = bspm_level2_regression(cons, regressors, varargin)
% BSPM_LEVEL2_REGRESSION
%
%   ARGUMENTS:
%       cons: contrast images from level 1
%       regressors: specified as follows: 
%           regressors(n).name = 'Regressor Name';
%           regressors(n).values = [1...n]; 
%

% ------------------------------------- Copyright (C) 2014 -------------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
def =   { ...
        'outdir',       [],                     ...
        'tag',          [],                     ...
        'implicitmask',  1,                     ...
        'mask',         '',                     ...
        'pctgroup',     [],                     ...
        'viewit',        0,                     ...
        'regidx',       [],                     ...
        'negativecon',   0,                     ...
        'nan2zero',      0,                     ...
        };
vals = setargs(def, varargin);
if nargin==0, mfile_showhelp; fprintf('\t= DEFAULT SETTINGS =\n'); disp(vals); return; end
fprintf('\n\t= CURRENT SETTINGS =\n'); disp(vals); 
if ischar(cons), cons1 = cellstr(cons1); end
if iscell(outdir), outdir = char(outdir); end
if iscell(mask), mask = char(mask); end
if ~isempty(mask)
    [~,mtag] = fileparts(mask); 
else
    mtag = 'noexplicitmask'; 
end

% | Regressors
nreg = length(regressors);

% | OUTPUT DIRECTORY
if isempty(outdir)
    
    % | Contrast Name
    hdr     = spm_vol(cons{1});
    idx     = strfind(hdr.descrip,':');
    cname   = hdr.descrip(idx+1:end);
    cname   = strtrim(regexprep(cname,'- All Sessions',''));
    
    % | Regressors
    regname = strcat('_', {regressors.name});
    regname = regexprep(regname, ' ', '');
    regname = strcat(regname{:});
    
    % | Analysis Name
    [p, level1name]  = fileparts(fileparts(cons{1})); 
    gadir       = fullfile(parentpath(cons), '_groupstats_', level1name);
    if tag
        gasubdir    = fullfile(gadir, sprintf('REGRESS%s_%s_N%d_%s_%s', tag, regname, length(cons), mtag, bspm_timestamp(1)));
    else
        gasubdir    = fullfile(gadir, sprintf('REGRESS%s_N%d_%s_%s', regname, length(cons), mtag, bspm_timestamp(1)));
    end
    outdir      = fullfile(gasubdir, cname);

    % | Make Directories
    if ~isdir(gadir), mkdir(gadir); end
    if ~isdir(gasubdir), mkdir(gasubdir); end
    
end
if ~isdir(outdir), mkdir(outdir); end

% | NAN2ZERO (IF APPLICABLE)
if nan2zero, bspm_batch_imcalc(cons, '', 'nan2zero'); end

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

% | FACTORIAL DESIGN SPECIFICATION
matlabbatch{1}.spm.stats.factorial_design.dir{1}            = outdir;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = cellstr(cons);
for i = 1:nreg
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).c      = regressors(i).values;
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).cname  = regressors(i).name; 
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).iCC    = 1;
end
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(nreg+1).c      = ones(length(cons), 1); 
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(nreg+1).cname  = 'Intercept';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(nreg+1).iCC    = 5;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint           = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.im                = implicitmask;
matlabbatch{1}.spm.stats.factorial_design.masking.em{1}             = mask; 
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit            = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no    = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm           = 1;

% | ESTIMATE
matlabbatch{2}.spm.stats.fmri_est.spmmat{1}        = fullfile(outdir,'SPM.mat');
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% | CONTRASTS
matlabbatch{3}.spm.stats.con.spmmat{1}               = fullfile(outdir,'SPM.mat');
matlabbatch{3}.spm.stats.con.delete                  = 1;
if isempty(regidx), regidx = 1:nreg; end
cnt = 0; 
for i = 1:length(regidx)
    cnt = cnt + 1;
    tmp = zeros(1,nreg); 
    tmp(regidx(i)) = 1; 
    matlabbatch{3}.spm.stats.con.consess{cnt}.tcon.name    = [regressors(regidx(i)).name '_POSITIVE']; 
    matlabbatch{3}.spm.stats.con.consess{cnt}.tcon.weights = tmp;
    matlabbatch{3}.spm.stats.con.consess{cnt}.tcon.convec  = tmp;
    matlabbatch{3}.spm.stats.con.consess{cnt}.tcon.sessrep = 'none';
    if negativecon
        cnt = cnt + 1; 
        matlabbatch{3}.spm.stats.con.consess{cnt}.tcon.name    = [regressors(regidx(i)).name '_NEGATIVE']; 
        matlabbatch{3}.spm.stats.con.consess{cnt}.tcon.weights = tmp*-1;
        matlabbatch{3}.spm.stats.con.consess{cnt}.tcon.convec  = tmp*-1;
        matlabbatch{3}.spm.stats.con.consess{cnt}.tcon.sessrep = 'none';
    end
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
