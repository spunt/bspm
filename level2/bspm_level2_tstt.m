function matlabbatch = bspm_level2_tstt(cons1, cons2, labels, outaffix, mask, outdir)
% BSPM_LEVEL2_TSTT
%
%   ARGUMENTS:
%       cons1: contrast images from level 1 for group 1
%       cons2: contrast images from level 1 for group 2
%       labels: labels for group 1 and group 2 (e.g., {'group1' 'group2'})
%       (optional) outaffix: string to affix to analysis directory name
%       (optional) mask: mask file to use
%

% -------------------------- Copyright (C) 2014 --------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014


% - check for existence of inputs and assign defaults - %
if nargin<6, outdir=[]; end
if nargin<5, mask=''; end
if nargin<4, outaffix='TSST'; end
if nargin<3, disp('USAGE: bspm_level2_tstt(cons1, cons2, labels, outaffix, outdir, mask)'); return; end

% - check formatting of inputs  - %
if ischar(cons1), cons1 = cellstr(cons1); end
if ischar(cons2), cons2 = cellstr(cons2); end
if iscell(outdir), outdir = char(outdir); end
if iscell(mask), mask = char(mask); end

% ------------------------------------------------------------

% - define output directory
if isempty(outdir)
    if iscell(outaffix), outaffix = char(outaffix); end
    [p n e] = fileparts(cons1{1});
    idx = strfind(p,'/');
    sname = p(1:idx(end-2)-1);
    aname = p(idx(end)+1:length(p));
    hdr = spm_vol(cons1{1});
    hdr = hdr.descrip;
    idx = strfind(hdr,':');
    cname = hdr(idx+1:end);
    cname = regexprep(cname,'- All Sessions','');
    cname = strtrim(cname);
    gname = [labels{1} '_-_' labels{2}];
    gadir = fullfile(sname,'_groupstats_',aname);
    outdir = fullfile(sname,'_groupstats_',aname,[gname '_' cname '_' outaffix]);
    if ~isdir(gadir)
        mkdir(gadir);
    end
    mkdir(outdir);
else
    if ~isdir(outdir), mkdir(outdir);end
end

% - fix end of image filename cell array
for i = 1:length(cons1), cons1(i) = cellstr([cons1{i} ',1']); end
for i = 1:length(cons2), cons2(i) = cellstr([cons2{i} ',1']); end

% - design
matlabbatch{1}.spm.stats.factorial_design.dir{1} = outdir;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = cellstr(cons1);
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = cellstr(cons2);
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
if ~isempty(mask)
    matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = strcat(mask,',1');
end
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
% - estimate
matlabbatch{2}.spm.stats.fmri_est.spmmat{1} = fullfile(outdir,'SPM.mat');
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
% - contrasts
matlabbatch{3}.spm.stats.con.spmmat{1} = fullfile(outdir,'SPM.mat');
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = [labels{1} '_-_' labels{2}];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = [labels{2} '_-_' labels{1}];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

% - run job
if nargout==0,  spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end


end
