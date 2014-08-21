function [] = bspm_Fcontrasts(analysis_dir, weights, delete_tag)
% BSPM_CONTRASTS
%
%   ARGUMENTS:
%       analysis_dir: directory containing SPM.mat
%       weights: contrast weights
%       delete_tag: tag for deleting existing contrasts: 0 for no (default), 1 for yes
%

% -------------------------------- Copyright (C) 2014 --------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<3
    delete_tag = 0;
elseif nargin<2
   disp('USAGE: bspm_contrasts(analysis_dir, weights, delete_tag)');
   return
end

name = 'Omnibus';

if ischar(analysis_dir)
    cellstr(analysis_dir);
end

for s = 1:length(analysis_dir)

spmmat = [analysis_dir{s} filesep 'SPM.mat'];
tmp = load(spmmat);

% build job
matlabbatch{1}.spm.stats.con.spmmat{1} = spmmat;
matlabbatch{1}.spm.stats.con.delete = delete_tag;   
matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = name;
matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'repl'; 
for r = 1:size(weights,1)
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.convec{r} = weights(r,:);
end

% run job
spm('defaults','fmri'); spm_jobman('initcfg');         
spm_jobman('run',matlabbatch);

end

end

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
