function bspm_level12(input)
% BSPM_LEVEL12
%
%   ARGUMENTS:
%
%       images:   functional images
%
%       general_info:   a structure with the following fields
%           analysis:   path for the analysis
%           TR:   repetition time (in seconds)
%           mt_res:   microtime resolution (number of time bins per scan)
%           mt_onset:   microtime onset
%           hpf:   high-pass filter cutoff to use (in seconds)
%           hrf_derivs:   temporal(1) and/or dispersion(2), e.g. [1 1] = use both
%           autocorrelation:   0=None, 1=AR(1), 2=Weighted Least Squares (WLS)
%           nuisance_file:   txt file with nuisance regressors (leave empty for none)
%           brainmask:   brainmask to use (leave empty for none)
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
%       contrasts:  a structure with the following fields
%           type:   'T' or 'F'
%           name:   string naming the contrast
%           weights:    vector of contrast weights
%           repl_tag:   tag to replicate weights across sessions (default=1)
%

% --------------------------------- Copyright (C) 2014 ---------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 1, error('No input!'); end
fn = {'images' 'general_info' 'runs' 'contrasts'};
[status, msg] = checkfields(input, fn);
if ~status, error(msg); end
images = input.images;
general_info = input.general_info;
runs = input.runs;
contrasts = input.contrasts;
if ~isfield(general_info, 'mt_res'), general_info.mt_res = 16; end
if ~isfield(general_info, 'mt_onset'), general_info.mt_onset = 1; end
if ~isfield(general_info, 'hrf_derivs'), general_info.hrf_derivs = [0 0]; end

% Session Non-Specific Parameters
mkdir(general_info.analysis);
matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.dir{1} = general_info.analysis;
matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.timing.units = 'secs';
matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.timing.RT = general_info.TR;
matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.timing.fmri_t = general_info.mt_res;
matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.timing.fmri_t0 = general_info.mt_onset; 
matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.bases.hrf.derivs(1) = general_info.hrf_derivs(1); % time derivative (0=no, 1=yes)
matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.bases.hrf.derivs(2) = general_info.hrf_derivs(2); % dispersion derivative (0=no, 1=yes)
matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.volt = 1;
matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.global = 'None';
matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.mask{1} = general_info.brainmask;
if general_info.autocorrelation==1
    matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.cvi = 'AR(1)';
elseif general_info.autocorrelation==2
    matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.cvi = 'wls';
else
    matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.cvi = 'none';
end

% Session Specific Parameters
nruns = length(runs);
if nruns==0
    runs.conditions = [];
    nruns = 1;
end
    
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
    matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(r).scans = cimages;
    for c = 1:length(conditions)
        matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(r).cond(c).name = conditions(c).name;
        matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(r).cond(c).onset = conditions(c).onsets;
        matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(r).cond(c).duration = conditions(c).durations;
        matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(r).cond(c).tmod = 0;
        if isfield(conditions(c), 'parameters')
            parameters = conditions(c).parameters;
            for p = 1:length(parameters)
                matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(r).cond(c).pmod(p).name = parameters(p).name;
                matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(r).cond(c).pmod(p).param = parameters(p).values;
                matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(r).cond(c).pmod(p).poly = 1;
            end
        end  
    end
    
    matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(r).multi{1} = '';
    matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(r).multi_reg{1} = nuisance;
    matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(r).hpf = general_info.hpf;

end
   
% Estimation Job
matlabbatch{2}.spm.tools.rwls.fmri_rwls_est.spmmat{1} = [general_info.analysis filesep 'SPM.mat'];         
matlabbatch{2}.spm.tools.rwls.fmri_rwls_est.method.Classical = 1;     

% Contrast Job
docon = 1;
if docon
    matlabbatch{3}.spm.stats.con.spmmat{1} = [general_info.analysis filesep 'SPM.mat'];
    matlabbatch{3}.spm.stats.con.delete = 1;   
    for c = 1:length(contrasts)
        
        if ~isfield(contrasts(c), 'repl_tag'), repl_tag = 1;
        else repl_tag = contrasts(c).repl_tag; end
        if repl_tag, repl_choice = 'repl'; else repl_choice = 'none'; end

        if strcmp(contrasts(c).type, 'T')
            matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = contrasts(c).name;
            matlabbatch{3}.spm.stats.con.consess{c}.tcon.convec = contrasts(c).weights;
            matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = repl_choice;
        else
            matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = contrasts(c).name;
            matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = repl_choice;
            weights = contrasts(c).weights;
            for r = 1:size(weights,1)
                matlabbatch{3}.spm.stats.con.consess{c}.fcon.convec{r} = weights(r,:);
            end
        end

    end
end

% run job
spm('defaults','fmri');
spm_jobman('run',matlabbatch);

end

 
 
 
 
