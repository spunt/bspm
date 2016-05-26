function matlabbatch = bspm_level2_ff_oneway(condirpat, factors, varargin)
% BSPM_LEVEL2_FF_ONEWAY
%
%   ARGUMENTS:
%
%       condirpat: pattern for finding level 1 dirs containing con images
%
%       factors: specified as follows: 
%           factors(n).name = 'Factor Name';
%           factors(n).idx2con = indices for con files corresponding to levels
%factors(1).idx2con = [1 2 3; 4 5 6]
% factors(1).name = 'WhyHow' 
%       covariates: specified as follows: 
%           covariates(n).name = 'Covariate Name';
% job = bspm_level2_ff_oneway(condirpat, factors, 'pctgroup', pctgroup, 'conweights', conweights);
%           covariates(n).values = [1...n]; 
%

% -------------------------- Copyright (C) 2014 --------------------------
%   Author: Bob Spunt
%   Affilitation: Caltech
%   Email: spunt@caltech.edu
%
%   $Revision Date: Aug_20_2014
def =   { ...
        'covariates',   [],                     ...
        'conpat',       'con*nii',              ...
        'conweights',   'auto',                 ...
        'implicitmask',  0,                     ...
        'mask',         '',                     ...
        'nan2zero',      1,                     ...
        'negativecon',   0,                     ...
        'outdir',       [],                     ...
        'pctgroup',     90,                     ...
        'tag',          [],                     ...
        'viewit',        0,                     ...
        };
vals = setargs(def, varargin);
if nargin<2, mfile_showhelp; fprintf('\t= DEFAULT SETTINGS =\n'); disp(vals); return; end
fprintf('\n\t= CURRENT SETTINGS =\n'); disp(vals);
if iscell(condirpat), cons1 = cellstr(condirpat); end
if ~regexp(condirpat, pwd), condirpat = fullfile(pwd, condirpat); end 
if iscell(outdir), outdir = char(outdir); end
if iscell(mask), mask = char(mask); end
if ~isempty(mask)
    [~,mtag] = fileparts(mask);
    mtag = upper(mtag);
else
    mtag = 'NOMASK'; 
end

% | GET CONS
conidx = sort(factors(1).idx2con(:));
[cpath, conlist] = grabfiles(condirpat, conpat, conidx, 1); 
if isempty(cpath), error('No contrasts found!'); end
[ucon, idx2ex, idx2sub] = unique(conlist);
conname = bspm_con2name(cpath(idx2ex));

% | BUILD MATRIX
conlistidx = find(ismember(conlist, ucon(idx2ex)));
cons    = cpath(conlistidx);
ncon    = length(cons);
ncell   = length(idx2ex);
nfact   = length(factors);
nsub    = length(cons)/ncell;
I       = ones(ncon, 4);
I(:,2)  = reshape(repmat(1:nsub,ncell,1), ncon, 1);
I(:,3)  = repmat(1:ncell, 1, nsub)';

% | NAN2ZERO (IF APPLICABLE)
if nan2zero, bspm_batch_imcalc(cons, '', 'nan2zero'); end

% | OUTPUT DIRECTORY
if isempty(outdir)
    
    % | Analysis Name
    [p, level1name] = fileparts(fileparts(cons{1})); 
    gadir           = fullfile(parentpath(cons), '_groupstats_', level1name);
    if isempty(pctgroup), pctgrouptag = 100; else pctgrouptag = pctgroup; end
    if tag
        gasubdir    = fullfile(gadir, sprintf('FF_%s_%s_N%d_PCTIN%d_%s_%s', tag, factors.name, nsub, pctgrouptag, mtag, bspm_timestamp(1)));
    else
        gasubdir    = fullfile(gadir, sprintf('FF_%s_N%d_PCTIN%d_%s_%s', factors.name, nsub, pctgrouptag, mtag, bspm_timestamp(1)));
    end
    outdir          = fullfile(gasubdir);

    % | Make Directories
    if ~isdir(gadir), mkdir(gadir); end
    
end
if ~isdir(outdir), mkdir(outdir); end

% | PCTGROUP (IF APPLICABLE)
if ~isempty(pctgroup)
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
cons = strcat(cons, ',1'); 
if ~isempty(mask), mask = strcat(mask, ',1'); end

% | FACTORIAL DESIGN SPECIFICATION
matlabbatch{1}.spm.stats.factorial_design.dir{1} = outdir;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'Subject';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
for i = 1:nfact
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(i+1).name     = factors(i).name;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(i+1).dept     = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(i+1).variance = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(i+1).gmsca    = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(i+1).ancova   = 0;
end
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.scans    = cellstr(cons);
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.imatrix  = I;

matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 2;
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.inter.fnums = [2 3];

% | MASKING & GLOBAL CALCULATIONS
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

if strcmpi(conweights, 'auto'), conweights = bspm_conweights(nlevel); end
if negativecon, conweights = [conweights; conweights*-1]; end
if ~isempty(conweights)
    matlabbatch{3}.spm.stats.con.spmmat{1}               = fullfile(outdir,'SPM.mat');
    matlabbatch{3}.spm.stats.con.delete                  = 1;
    for c = 1:size(conweights,1);
        matlabbatch{3}.spm.stats.con.consess{c}.tcon.name    = bspm_conweights2names(conweights(c,:), conname);
        matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = conweights(c,:);
        matlabbatch{3}.spm.stats.con.consess{c}.tcon.convec  = conweights(c,:);
        matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
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
