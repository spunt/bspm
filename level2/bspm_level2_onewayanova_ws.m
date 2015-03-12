function matlabbatch = bspm_level2_onewayanova_ws(cellcons, outname, mask)
% BSPM_LEVEL2_ONEWAYANOVA_WS
%
%   ARGUMENTS:
%       cons: contrast images from level 1
%       outname: name for directory in _groupstats_/<analysis_name>
%       (optional) mask: mask file to use (default = none)
%

% ----------------------- Copyright (C) 2014 -----------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<2, disp('USAGE: bspm_level2_onewayanova_ws(cons, outname, mask)'); return
elseif nargin<3, mask = ''; end
if iscell(mask), mask = char(mask); end
if iscell(outname), outname = char(outname); end

% define output directory
tmp = cellcons{1}; [p n e] = fileparts(tmp{1});
idx = strfind(p,'/');
sname = p(1:idx(end-2)-1);
aname = p(idx(end)+1:length(p));
gadir = fullfile(sname,'_groupstats_',aname);
if ~isdir(gadir), mkdir(gadir); end
outdir = fullfile(sname,'_groupstats_',aname,outname);
mkdir(outdir);

% re-arrange by subject
ncond = length(cellcons);
nsub = length(cellcons{1});
for s = 1:nsub
    for c = 1:ncond
        subcon{s}(c) = cellcons{c}(s);
    end
end

% get names
for c = 1:ncond
    hdr = spm_vol(subcon{1}{c});
    hdr = hdr.descrip;
    idx = strfind(hdr,':');
    cname = hdr(idx+1:end);
    cname = regexprep(cname,'- All Sessions','');
    connam{c} = strtrim(cname);
end

% build job
for s = 1:length(subcon)
    cons = subcon{s};
    for i = 1:length(cons), cons(i) = cellstr([cons{i} ',1']); end
    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(s).scans = cellstr(cons);
    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(s).conds = 1:ncond;
end
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.dept = 1; % independence: 0 = Yes, 1 = No
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.variance = 1; % variance: 0 = Equal, 1 = Unequal
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1; % implicit masking
if ~isempty(mask)
    matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = strcat(mask,',1'); % explicit mask
end
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{1}.spm.stats.factorial_design.dir{1} = outdir;
matlabbatch{2}.spm.stats.fmri_est.spmmat{1} = fullfile(outdir,'SPM.mat');
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat{1} = fullfile(outdir,'SPM.mat');
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = [connam{1} '_-_' connam{2}];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = [connam{2} '_-_' connam{1}];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

% run job
if nargout==0,  spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end


end
