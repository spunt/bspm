function [] = bspm_level1_specify_design(images, general_info, runs)
% BSPM_LEVEL1_SPECIFY_DESIGN
%
%   USAGE: bspm_level1_specify_design(images, general_info, runs)
%
%   ARGUMENTS:
%
%       images:   functional images
%
%       general_info:   a structure with the following fields
%           analysis:   path for the analysis
%           TR:   repetition time (in seconds)
%           hpf:   high-pass filter cutoff to use (in seconds)
%           autocorrelation:   0=None, 1=AR(1)
%           nuisance_file:   txt file with nuisance regressors (leave empty for none)
%
%       runs:   a structure with the following fields
%           conditions:
%               name:   string naming the condition
%               onsets:   onsets
%               durations:   durations
%               parameters:   for building parametric modulators (leave empty for none)
%                   name:   string naming the paramter
%                   values:   parameter values (assumed to be orthogonalized)
%

% --------------------------------- Copyright (C) 2014 ---------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<3, disp('bspm_level1_specify_design(images, general_info, runs)'); return; end

% Session Non-Specific Parameters
mkdir(general_info.analysis);
matlabbatch{1}.spm.stats.fmri_design.dir{1} = general_info.analysis;
matlabbatch{1}.spm.stats.fmri_design.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_design.timing.RT = general_info.TR;
matlabbatch{1}.spm.stats.fmri_design.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_design.timing.fmri_t0 = 1; 
matlabbatch{1}.spm.stats.fmri_design.bases.hrf.derivs(1) = 0; % time derivative (0=no, 1=yes)
matlabbatch{1}.spm.stats.fmri_design.bases.hrf.derivs(2) = 0; % dispersion derivative (0=no, 1=yes)
matlabbatch{1}.spm.stats.fmri_design.volt = 1;
matlabbatch{1}.spm.stats.fmri_design.global = 'None';
if general_info.autocorrelation==1
    matlabbatch{1}.spm.stats.fmri_design.cvi = 'AR(1)';
else
    matlabbatch{1}.spm.stats.fmri_design.cvi = 'none';
end

% Session Specific Parameters
nruns = length(runs);
for r = 1:nruns
    if nruns==1
        cimages = images;
        conditions = runs.conditions;
        nuisance = general_info.nuisance_file;
    else
        cimages = images{r};
        conditions = runs(r).conditions;
        nuisance = general_info.nuisance_file{r};
    end
    matlabbatch{1}.spm.stats.fmri_design.sess(r).nscan = length(cimages);
    for c = 1:length(conditions)
        matlabbatch{1}.spm.stats.fmri_design.sess(r).cond(c).name = conditions(c).name;
        matlabbatch{1}.spm.stats.fmri_design.sess(r).cond(c).onset = conditions(c).onsets;
        matlabbatch{1}.spm.stats.fmri_design.sess(r).cond(c).duration = conditions(c).durations;
        matlabbatch{1}.spm.stats.fmri_design.sess(r).cond(c).tmod = 0;
        if isfield(conditions(c), 'parameters')
            parameters = conditions(c).parameters;
            for p = 1:length(parameters)
                matlabbatch{1}.spm.stats.fmri_design.sess(r).cond(c).pmod(p).name = parameters(p).name;
                matlabbatch{1}.spm.stats.fmri_design.sess(r).cond(c).pmod(p).param = parameters(p).values;
                matlabbatch{1}.spm.stats.fmri_design.sess(r).cond(c).pmod(p).poly = 1;
            end
        end  
    end
    matlabbatch{1}.spm.stats.fmri_design.sess(r).multi{1} = '';
    matlabbatch{1}.spm.stats.fmri_design.sess(r).multi_reg{1} = nuisance;
    matlabbatch{1}.spm.stats.fmri_design.sess(r).hpf = general_info.hpf;
end
   
% run job
spm('defaults','fmri'); spm_jobman('initcfg');         
spm_jobman('run',matlabbatch);

end

 
 
 
 
