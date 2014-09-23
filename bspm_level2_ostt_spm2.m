function [] = bspm_level2_ostt_spm2(input)
% BSPM_LEVEL2_OSTT
%
%   ARGUMENTS:
%       cons: contrast images from level 1
%       mask: mask file to use (default = none)
%       implicitTAG: 0 = no implicit masking; 1 = yes (default)
%

% ---------------------------------- Copyright (C) 2014 ----------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin < 1, error('No input!'); end
if ~isfield(input, 'mask'), input.mask = ''; end
if ~isfield(input, 'implicitTAG'), input.implicitTAG = 1; end
if ~isfield(input, 'outdir'), input.outdir = []; end
fn = {'cons' 'mask' 'implicitTAG' 'outdir'};
[status, msg] = checkfields(input, fn);
if ~status, error(msg); end

cons = input.cons; 
mask = input.mask; 
implicitTAG = input.implicitTAG; 
outdir = input.outdir; 
if ischar(cons), cons = cellstr(cons); end
if iscell(mask), mask = char(mask); end
if ~isempty(mask)
    [tmpp,maskstr] = fileparts(mask); 
else
    maskstr = 'NoMask'; 
end
name = sprintf('N=%d_%s', length(cons), maskstr); 

% define output directory
if isempty(outdir)
    
    if iscell(name), name = char(name); end
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
    gasubdir = fullfile(gadir,name);
    outdir = fullfile(gasubdir,cname);
    if ~isdir(gadir), mkdir(gadir); end
    if ~isdir(gasubdir), mkdir(gasubdir); end
    mkdir(outdir);
    
else
    if ~isdir(outdir), mkdir(outdir);end
end

% fix end of image filename cell array
for i = 1:length(cons)
    cons(i) = cellstr([cons{i} ',1']);
end

% build job
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = cellstr(cons);
matlabbatch{1}.spm.stats.factorial_design.masking.im = implicitTAG;
if ~isempty(mask)
    matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = strcat(mask,',1');
end
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{1}.spm.stats.factorial_design.dir{1} = outdir;
matlabbatch{2}.spm.stats.fmri_est.spmmat{1} = fullfile(outdir,'SPM.mat');
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat{1} = fullfile(outdir,'SPM.mat');
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Positive';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Negative';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = -1;
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

% run job
spm('defaults','fmri'); spm_jobman('initcfg');         
spm_jobman('run',matlabbatch);

end

 
 
 
 
