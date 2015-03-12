function matlabbatch = 
% BSPM_LEVEL2_REGRESSION
%
%   ARGUMENTS:
%       cons: contrast images from level 1
%       regressor: vector of values
%       regname: name for regressor
%       outprefix: prefix for directory  in _groupstats_/<analysis_name>/<regname>
%       (optional) mask: mask file to use (default = none)
%       (optional) implicitTAG: 0 = no implicit masking; 1 = yes (default)
%

% ------------------------------------- Copyright (C) 2014 -------------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<6, implicitTAG=1; end
if nargin<5, mask = ''; end
if nargin<4, disp('USAGE: bspm_level2_regression(cons, regressor, regname, outprefix, mask, implicitTAG)'); return; end

% make sure image names
if ischar(cons)
    cons = cellstr(cons);
end
if iscell(mask)
    mask = char(mask);
end
if iscell(outprefix)
    outprefix = char(outprefix);
end

% define output directory
[p n e] = fileparts(cons{1});
idx = strfind(p,'/');
sname = p(1:idx(end-2)-1);
aname = p(idx(end)+1:length(p));
hdr = spm_vol(cons{1});
hdr = hdr.descrip;
idx = strfind(hdr,':');
cname = hdr(idx+1:end);
cname = regexprep(cname,'- All Sessions','');
cname = strtrim(cname);
gadir = fullfile(sname,'_groupstats_',aname);
regdir = fullfile(gadir, regname);
outdir = fullfile(regdir,[outprefix '_' cname]);
if ~isdir(gadir)
    mkdir(gadir);
end
if ~isdir(regdir)
    mkdir(regdir);
end
mkdir(outdir);

% fix end of image filename cell array
for i = 1:length(cons)
    cons(i) = cellstr([cons{i} ',1']);
end

% build job
matlabbatch{1}.spm.stats.factorial_design.dir{1} = outdir;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = cellstr(cons);
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = regressor;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = regname;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {}); %
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = implicitTAG;
if ~isempty(mask)
    matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = strcat(mask,',1');
end
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat{1} = fullfile(outdir,'SPM.mat');
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat{1} = fullfile(outdir,'SPM.mat');
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Positive';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [0 1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Negative';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [0 -1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

% make outdir
mkdir(outdir)

% run job
if nargout==0,  spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end


end
