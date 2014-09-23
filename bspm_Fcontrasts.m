function [] = bspm_Fcontrasts(analysis_dirs, weights, delete_tag)
% BSPM_FCONTRASTS
%
%   USAGE:  bspm_Fcontrasts(analysis_dirs, weights, [delete_tag])
%
%   ARGUMENTS:
%      analysis_dirs: directories containing SPM.mat
%      weights: contrast weights
%      delete_tag: delete existing contrasts: 0 for no (default), 1 for yes
%

% ---------------------------- Copyright (C) 2014 -------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Sep_23_2014
if nargin<2, error('USAGE: bspm_Fcontrasts(analysis_dirs, weights, delete_tag)'); end
if nargin<3, delete_tag = 0; end
if ischar(analysis_dirs), cellstr(analysis_dirs); end
for s = 1:length(analysis_dirs)

    spmmat = [analysis_dirs{s} filesep 'SPM.mat'];
    tmp = load(spmmat);

    % build job
    matlabbatch{1}.spm.stats.con.spmmat{1} = spmmat;
    matlabbatch{1}.spm.stats.con.delete = delete_tag;   
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'Omnibus';
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'repl'; 
    for r = 1:size(weights,1)
        matlabbatch{1}.spm.stats.con.consess{1}.fcon.convec{r} = weights(r,:);
    end

    % run job
    spm('defaults','fmri'); spm_jobman('initcfg');         
    spm_jobman('run',matlabbatch);

end
end

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
