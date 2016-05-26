function matlabbatch = bspm_level1(images, general_info, runs, contrasts)
% BSPM_LEVEL1
%
%   USAGE: matlabbatch = bspm_level1(images, general_info, runs, contrasts)
%
%   ARGUMENTS:
%
%       images:                 nruns x 1 cell array with functional image files
% 
%       general_info - structure with the following fields:
%           analysis:           path for the analysis
%           is4D:               0 = 3D images (default), 1 = 4D image
%           TR:                 repetition time (in seconds)
%           mt_res:             microtime resolution (number of time bins per scan)
%           mt_onset:           microtime onset
%           hpf:                high-pass filter cutoff to use (in seconds)
%           hrf_derivs:         temporal(1) and/or dispersion(2), e.g. [1 1] = use both
%           autocorrelation:    0=None, 1=AR(1), 2=Weighted Least Squares (WLS), 3=FAST
%           nuisance_file:      txt file with nuisance regressors (leave empty for none)
%           brainmask:          brainmask to use (leave empty for none)
%           maskthresh:         masking threshold (proportion of globals), default = -Inf
%           orth:               0 = do not orthogonalise pmod regressors (default), 1 = do
% 
%       runs - structure with the following fields
%           conditions:
%               name:           string naming the condition
%               onsets:         onsets
%               durations:      durations
%               parameters:     for within-condition parametric modulators (leave empty for none)
%                   name:       string naming the paramter
%                   values:     parameter values
%           floatingpm          for between-condition parametric modulators
%               name:       string naming the parameter
%               onsets:     onsets of modulated events
%               durations:  durations of modulated events
%               values:     parameter values
%           regressors:         for regressors not convolved with the HRF (e.g., motion)
%               name:           string naming the regressors
%               values:         regressors values (1 x nscans)
% 
%       contrasts - structure with the following fields
%           type:               'T' or 'F'
%           name:               string naming the contrast
%           weights:            vector of contrast weights
%           repl_tag:           tag to replicate weights across sessions (default = 1)
%

% --------------------------------- Copyright (C) 2014 ---------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
docon = 1; 
if nargin==1 
    input         = images;
    fn            = {'images' 'general_info' 'runs'};
    [status, msg] = checkfields(input, fn);
    if ~status, error(msg); end
    images        = input.images;
    general_info  = input.general_info;
    runs          = input.runs;
    if checkfields(input, 'contrasts'); 
        contrasts = input.contrasts;
    else
        docon     = 0; 
    end
else
    if nargin < 3, mfile_showhelp; return; end
    if nargin < 4, docon = 0; end
end

% | Assign Defaults
% | =======================================================================
if ~isfield(general_info, 'is4D'), general_info.is4D = 0; end
if ~isfield(general_info, 'mt_res'), general_info.mt_res = 16; end
if ~isfield(general_info, 'mt_onset'), general_info.mt_onset = 8; end
if ~isfield(general_info, 'hrf_derivs'), general_info.hrf_derivs = [0 0]; end
if ~isfield(general_info, 'maskthresh'), general_info.maskthresh = 0.8; end
if ~isfield(general_info, 'orth'), general_info.orth = 0; end
if ~isfield(general_info, 'write_residuals'), general_info.write_residuals = 0; end
if isempty(general_info.brainmask), general_info.brainmask = ''; end
if isempty(general_info.nuisance_file), general_info.nuisance_file = ''; end
if ~exist(general_info.analysis, 'dir'), mkdir(general_info.analysis); end

% | Specify Design
% | =======================================================================

% | Session Non-Specific Parameters
spec.dir{1}                 = general_info.analysis;
spec.timing.units           = 'secs';
spec.timing.RT              = general_info.TR;
spec.timing.fmri_t          = general_info.mt_res;
spec.timing.fmri_t0         = general_info.mt_onset;
spec.fact                   = struct('name', {}, 'levels', {});
spec.bases.hrf.derivs(1)    = general_info.hrf_derivs(1); % time derivative (0=no, 1=yes)
spec.bases.hrf.derivs(2)    = general_info.hrf_derivs(2); % dispersion derivative (0=no, 1=yes)
spec.volt                   = 1;
spec.global                 = 'None';
spec.mask{1}                = general_info.brainmask;
spec.mthresh                = general_info.maskthresh; 
autocorropt                 = {'none' 'AR(1)' 'wls' 'FAST'};
spec.cvi                    = autocorropt{general_info.autocorrelation+1}; 

% | Session Specific Parameters
nruns = length(runs);
if nruns==0, runs.conditions = []; nruns = 1; end
for r = 1:nruns
    regressors = [];
    floatingpm = []; 
    if nruns==1
        cimages = images;
        conditions = runs.conditions;
        nuisance = general_info.nuisance_file;
        if isfield(runs, 'regressors'), regressors = runs.regressors; end
        if isfield(runs, 'floatingpm'), floatingpm = runs.floatingpm; end
    else
        cimages = images{r};
        conditions = runs(r).conditions;
        nuisance = general_info.nuisance_file{r};
        if isfield(runs(r), 'regressors'), regressors = runs(r).regressors; end
    end
    if general_info.is4D
        cimages = bspm_expand4D(cimages);
    end
    spec.sess(r).scans = cimages;
    for c = 1:length(conditions)
        spec.sess(r).cond(c).name     = conditions(c).name;
        spec.sess(r).cond(c).onset    = conditions(c).onsets;
        spec.sess(r).cond(c).duration = conditions(c).durations;
        spec.sess(r).cond(c).tmod     = 0;
        if isfield(conditions(c), 'parameters')
            parameters = conditions(c).parameters;
            for p = 1:length(parameters)
                spec.sess(r).cond(c).pmod(p).name  = parameters(p).name;
                spec.sess(r).cond(c).pmod(p).param = parameters(p).values;
                spec.sess(r).cond(c).pmod(p).poly  = 1;
            end
            spec.sess(r).cond(c).orth = general_info.orth;
        end  
    end
    rc = 0; 
    if ~isempty(floatingpm)
        for p = 1:length(floatingpm)
            x = bspm_make_regressor(length(cimages), general_info.TR, floatingpm(p).onsets, ... 
                floatingpm(p).durations, 'PMods', floatingpm(p).values, 'TRbin', general_info.mt_res, ...
                'TRons', general_info.mt_onset, 'Derivatives', sum(general_info.hrf_derivs));
            rc = rc + 1; 
            spec.sess(r).regress(rc).name   = floatingpm(p).name;
            spec.sess(r).regress(rc).val    = x(:,2); 
            if general_info.hrf_derivs(1)
                rc = rc + 1; 
                spec.sess(r).regress(rc).name   = strcat([floatingpm(p).name, '_xTD']); 
                spec.sess(r).regress(rc).val    = x(:,3); 
            end
        end
    end
    if ~isempty(regressors)
        for p = 1:length(regressors)
            rc = rc + 1; 
            spec.sess(r).regress(rc).name = regressors(p).name;
            spec.sess(r).regress(rc).val  = regressors(p).values;
        end
    end  
    spec.sess(r).multi{1}     = '';
    spec.sess(r).multi_reg{1} = nuisance;
    spec.sess(r).hpf          = general_info.hpf;
end
    
% | Estimation Job
% | =======================================================================
est.spmmat{1}        = [general_info.analysis filesep 'SPM.mat'];
est.write_residuals  = general_info.write_residuals;
est.method.Classical = 1;

% | Check for Use of RobustWLS Toolbox
% | =======================================================================
if strcmpi(spec.cvi, 'wls')
    matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec = spec; 
    matlabbatch{2}.spm.tools.rwls.fmri_rwls_est = est;
else
    matlabbatch{1}.spm.stats.fmri_spec = spec; 
    matlabbatch{2}.spm.stats.fmri_est = est;
end

% | Contrast Job
% | =======================================================================
if docon
    matlabbatch{3}.spm.stats.con.spmmat{1}  = fullfile(general_info.analysis, 'SPM.mat');
    matlabbatch{3}.spm.stats.con.delete     = 1;   
    for c = 1:length(contrasts)
        if ~isfield(contrasts(c), 'repl_tag'), repl_tag = 1;
        else repl_tag = contrasts(c).repl_tag; end
        if repl_tag, repl_choice = 'repl'; else repl_choice = 'none'; end
        if strcmpi(contrasts(c).type, 'T')
            matlabbatch{3}.spm.stats.con.consess{c}.tcon.name       = contrasts(c).name;
            matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights    = contrasts(c).weights; % SPM12
            matlabbatch{3}.spm.stats.con.consess{c}.tcon.convec     = contrasts(c).weights; % SPM8
            matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep    = repl_choice;
        else
            matlabbatch{3}.spm.stats.con.consess{c}.fcon.name       = contrasts(c).name;
            matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights    = contrasts(c).weights; % SPM12
            matlabbatch{3}.spm.stats.con.consess{c}.fcon.convec{1}  = contrasts(c).weights; % SPM8
            matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep    = repl_choice;
        end
    end
end

% | Run job (only if no output arguments requested)
% | =======================================================================
if nargout==0, spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end

end
 
 
 
 
