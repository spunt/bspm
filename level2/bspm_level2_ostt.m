function matlabbatch = bspm_level2_ostt(cons, varargin)
% BSPM_LEVEL2_OSTT
%
%   USAGE: matlabbatch = bspm_level2_ostt(cons, varargin)
%
%         cons: contrast images from level 1
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

% ---------------------------------- Copyright (C) 2014 ----------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
def = { 'outdir',       [],     ... 
        'implicit',     1,      ...
        'mask',         '',     ...
        'pctgroup',     [],     ...
        'viewit',       0,      ...
        'nan2zero',     0 };
vals = setargs(def, varargin);
if nargin==0, mfile_showhelp; fprintf('\t= DEFAULT SETTINGS =\n'); disp(vals); return; end
fprintf('\n\t= CURRENT SETTINGS =\n'); disp(vals); 


if ischar(cons), cons = cellstr(cons); end
if iscell(mask), mask = char(mask); end
if ~isempty(mask)
    [~,mtag] = fileparts(mask); 
else
    mtag = 'noexplicitmask'; 
end

% | OUTPUT DIRECTORY
if isempty(outdir)
    
    % | Contrast Name
    hdr = spm_vol(cons{1});
    idx = strfind(hdr.descrip,':');
    cname = hdr.descrip(idx+1:end);
    cname = strtrim(regexprep(cname,'- All Sessions',''));
    
    % | Analysis Name
    [p, level1name]  = fileparts(fileparts(cons{1})); 
    gadir       = fullfile(parentpath(cons), '_groupstats_', level1name); 
    gasubdir    = fullfile(gadir, sprintf('OSTT_N=%d_%s_%s', length(cons), mtag, bspm_timestamp(1)));
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
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans      = cellstr(cons);
matlabbatch{1}.spm.stats.factorial_design.masking.im        = implicit;
matlabbatch{1}.spm.stats.factorial_design.masking.em{1}     = mask; 

% | FACTORIAL DESIGN SPECIFICATION
matlabbatch{2}.spm.stats.fmri_est.spmmat{1} = fullfile(outdir,'SPM.mat');
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% | CONTRASTS
matlabbatch{3}.spm.stats.con.spmmat{1} = fullfile(outdir,'SPM.mat');
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Positive';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Negative';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1;
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

% | RUN IF NO OUTPUT ARGS SPECIFIED
if nargout==0
    bspm_runbatch(matlabbatch);
    if viewit
        cd(outdir); 
        bspmview; 
    end
end

end
