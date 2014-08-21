function [] = bspm_level2_onewayanova(cellcons, outname, mask)
% BSPM_LEVEL2_ONEWAYANOVA
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

if nargin<2, disp('USAGE: bspm_level2_ostt(cons, outprefix, mask)'); return
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

% build job
for cell = 1:length(cellcons)
    cons = cellcons{cell};
    % get condition name
    hdr = spm_vol(cons{1});
    hdr = hdr.descrip;
    idx = strfind(hdr,':');
    cname = hdr(idx+1:end);
    cname = regexprep(cname,'- All Sessions','');
    connam{cell} = strtrim(cname);
    for i = 1:length(cons), cons(i) = cellstr([cons{i} ',1']); end
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(cell).scans = cellstr(cons);
end
matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = 1; % independence: 0 = Yes, 1 = No
matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = 1; % variance: 0 = Equal, 1 = Unequal
matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1; % implicit masking
if ~isempty(mask)
    matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = strcat(mask,',1'); % explicit mask
end
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{1}.spm.stats.factorial_design.dir{1} = outdir;
matlabbatch{2}.spm.stats.fmri_est.spmmat{1} = fullfile(outdir,'SPM.mat');
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat{1} = fullfile(outdir,'SPM.mat');
basevec = zeros(1,length(connam));
for c = 1:length(connam)
    cvec = basevec; cvec(c) = 1;
    matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = connam{c};
    matlabbatch{3}.spm.stats.con.consess{c}.tcon.convec = cvec;
    matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
end
matlabbatch{3}.spm.stats.con.delete = 1;

% run job
spm('defaults','fmri'); spm_jobman('initcfg');         
spm_jobman('run',matlabbatch);

end



 
 
 
 
