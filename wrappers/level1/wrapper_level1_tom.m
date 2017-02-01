function matlabbatch = wrapper_level1_tom(covidx, varargin)
% matlabbatch = wrapper_level1_tom(covidx, varargin)
%
% To show default settings, run without any arguments.
%
%     COVIDX 
%       01 - Duration          
%       02 - No Response
% 

% | SET DEFAULTS AND PARSE VARARGIN
% | ===========================================================================
defaults = {
            'studydir',         '/Users/bobspunt/Drive/Writing/Empirical/ASD/data',            ...
            'HPF',              128,            ...
            'armethod',         1,              ... 
            'nuisancepat',      'bad*txt',      ...
            'epipat',           'sw*nii',       ...
            'subid',            'RA*',          ...
            'runid',            'EP*TOM*',      ...
            'behavid',          'tom_*mat',     ...
            'basename',         'TOM_s8w3',     ...
            'brainmask',        bspm_brainmask, ...
            'fcontrast',        1,              ...
            'nskip',            2,              ...
            'runtest',          0,              ...
            'TR',               1,              ...
            'is4D',             1               ...
             };
vals = setargs(defaults, varargin);
if nargin==0, mfile_showhelp; fprintf('\t= DEFAULT SETTINGS =\n'); disp(vals); return; end
fprintf('\n\t= CURRENT SETTINGS =\n'); disp(vals); 

% | PATHS
% | ===========================================================================
if strfind(pwd,'/home/spunt'), studydir = '/home/spunt/data/asd'; end
[subdir, subnam] = files([studydir filesep subid]);

% | ANALYSIS NAME
% | ===========================================================================
armethodlabels  = {'NoAR1' 'AR1' 'WLS'};
covnames        = {'Duration' 'NR'};
pmnames         = regexprep(covnames(covidx), '_', '');
pmstr           = sprintf(repmat('_%s', 1, length(pmnames)), pmnames{:}); pmstr(1)= [];
analysisname    = sprintf('%s_Pmodby_%s_%s_%ds_%s', basename, ...
                        pmstr, armethodlabels{armethod + 1}, HPF, bob_timestamp);
printmsg(analysisname, 'msgtitle', 'Analysis Name');

% | IMAGING PARAMETERS
% | ========================================================================
adjons          = TR*nskip;
                                    
% | RUNTIME OPTIONS
% | ===========================================================================           
if runtest, subdir = subdir(1); end

% | SUBJECT LOOP
% | ===========================================================================
matlabbatch = []; 
for s = 1:length(subdir)
    
    % | Check Subject and Define Folders
    % | ========================================================================
    rundir      = files([subdir{s} filesep 'raw' filesep runid]);
    if isempty(rundir), printmsg('Valid run directory not found, moving on...', 'msgtitle', subnam{s}); continue; end
    analysisdir = fullfile(subdir{s}, 'analysis', analysisname); 
    if any([exist(fullfile(analysisdir, 'mask.img'), 'file') exist(fullfile(analysisdir, 'mask.nii'), 'file')])
        printmsg('Level 1 job probably already estimated, moving on...', 'msgtitle', subnam{s}); continue; 
    end
    printmsg(sprintf('Building Level 1 Job for %d Runs', length(rundir)), 'msgtitle', subnam{s}); 

    % | Behavioral and Nuisance Regressor Files
    % | ========================================================================
    nuisance    = files([subdir{s} filesep 'raw' filesep runid filesep nuisancepat]);
    behav       = files([subdir{s} filesep 'behav' filesep behavid]);
    
    % | Run Loop
    % | ========================================================================
    images          = cell(size(rundir)); 
    for r = 1:length(rundir)

        % | Data for Current Run
        % | =====================================================================
        images{r}   = files([rundir{r} filesep epipat]);
        if isempty(images{r})
            error('\nImage data not found! Failed search pattern:\n%s', [rundir{r} filesep epipat]); 
        end
        load(behav{r});
        Seeker(:,6:7) = Seeker(:,6:7) - adjons;
        condlabels = {'Belief' 'Photo'};
        % | Columns for Seeker
        % | =====================================================================
        % 1 - TRIAL #
        % 2 - CONDITION (1=Belief, 2=Photo)
        % 5 - STIMULUS IDX
        % 6 - STORY ONSET
        % 7 - QUESTION ONSET
        % 8 - RESPONSE
        % 9 - RT TO QUESTION ONSET
        % 10 - BLOCK DURATION

        % | Conditions
        % | =====================================================================
        for c = 1:length(condlabels)
            runs(r).conditions(c).name      = condlabels{c};
            runs(r).conditions(c).onsets    = Seeker(Seeker(:,2)==c,6);
            runs(r).conditions(c).durations = Seeker(Seeker(:,2)==c,10);
        end
        
        % | Floating Parametric Modulators
        % | =====================================================================
        allpm           = [Seeker(:,10) double(Seeker(:,8)==0)]; 
        modelpm         = allpm(:,covidx);
        modelpmnames    = pmnames; 
        novaridx = find(nanstd(modelpm)==0);
        if ~isempty(novaridx), modelpm(:,novaridx) = []; modelpmnames(novaridx) = []; end
        for p = 1:length(modelpmnames)
            runs(r).floatingpm(p).name      = modelpmnames{p};
            runs(r).floatingpm(p).onsets    = Seeker(:,6); 
            runs(r).floatingpm(p).durations = Seeker(:,10);
            runs(r).floatingpm(p).values    = modelpm(:,p);
        end
        
    end
    if length(rundir)==1 
        images = images{1}; 
        if iscell(nuisance), nuisance = nuisance{1}; end
    end

    % | General Information
    % | ========================================================================
    general_info.analysis           = analysisdir; 
    general_info.is4D               = is4D; 
    general_info.TR                 = TR;
    general_info.hpf                = HPF;
    general_info.autocorrelation    = armethod;
    general_info.nuisance_file      = nuisance;
    general_info.brainmask          = brainmask;
    general_info.hrf_derivs         = [0 0];
    general_info.mt_res             = 16; 
    general_info.mt_onset           = 1;

    % | Contrasts
    % | ========================================================================
    ncond   = length(condlabels);
    w1      = eye(ncond);
    w2      = [1 -1];
    weights = [w1; w2]; 
    ncon    = size(weights,1);
    for c = 1:ncon
        contrasts(c).type       = 'T';
        contrasts(c).weights    = weights(c,:);
        contrasts(c).name       = bspm_conweights2names(weights(c,:), condlabels);
    end
    if fcontrast
        contrasts(ncon+1).type      = 'F';
        contrasts(ncon+1).name      = 'Omnibus';
        contrasts(ncon+1).weights   = eye(ncond);
    end

    % | Make Job
    % | ========================================================================
    matlabbatch = [matlabbatch bspm_level1(images, general_info, runs, contrasts)]; 

    % | Cleanup Workspace
    % | ========================================================================
    clear general_info runs contrasts b modelpm modelpmnames

end
end
% =========================================================================
% * SUBFUNCTIONS
% =========================================================================
function mfile_showhelp(varargin)
% MFILE_SHOWHELP
ST = dbstack('-completenames');
if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
eval(sprintf('help %s', ST(2).file));  
end

 
 
 
 
 
 
 
 
