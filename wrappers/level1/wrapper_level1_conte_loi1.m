function allinput = wrapper_level1_conte_loi1(varargin)
% matlabbatch = wrapper_level1_loi1(varargin)
%
% To show default settings, run without any arguments.
%
%

% | SET DEFAULTS AND PARSE VARARGIN
% | ===========================================================================
defaults = {
            'studydir',         '/Users/bobspunt/Documents/fmri/conte', ...
            'covidx',           [],             ...
            'epifname',         [],             ...
            'HPF',              100,            ...
            'armethod',         2,              ...
            'nuisancepat',      'bad*txt',      ...
            'epipat',           'swbua*nii',    ...
            'subid',            'CC*',          ...
            'runid',            'EP*LOI_1*',    ...
            'behavid',          'loi1*mat',     ...
            'basename',         'LOI1_s6w2',    ...
            'brainmask',        '',             ...
            'fcontrast',        1,              ...
            'nskip',            4,              ...
            'runtest',          0,              ...
            'is4D',             1,              ...
            'TR',               1,              ...
            'yesnokeys',        [1 2]           ...
             };
vals = setargs(defaults, varargin);
if nargin==0, mfile_showhelp; fprintf('\t= DEFAULT SETTINGS =\n'); disp(vals); return; end
fprintf('\n\t= CURRENT SETTINGS =\n'); disp(vals);

% | PATHS
% | ===========================================================================
if strfind(pwd,'/data/spunt'), studydir = '/data/spunt/Documents/fmri/conte'; end
[subdir, subnam] = files([studydir filesep subid]);

% | EPI FNAME
% | ===========================================================================
if ~isempty(epifname)
    epifname = fullfile(studydir, epifname);
    if exist(epifname, 'file')
        fnepi = load(epifname);
    else
        disp('epifname could not be found!');
        fnepi = [];
    end
else
    fnepi = [];
end

% | ANALYSIS NAME
% | ===========================================================================
armethodlabels  = {'NoAR1' 'AR1' 'WLS'};
covnames        = {'Duration' 'Errors' 'FoilErrors'};
if ~isempty(covidx)
    pmnames         = regexprep(covnames(covidx), '_', '');
    pmstr           = sprintf(repmat('_%s', 1, length(pmnames)), pmnames{:}); pmstr(1)= [];
else
    pmstr = 'None';
end
analysisname  = sprintf('%s_%s_%ds_%s', basename, ...
                armethodlabels{armethod + 1}, HPF, bob_timestamp);
printmsg(analysisname, 'msgtitle', 'Analysis Name');

% | IMAGING PARAMETERS
% | ========================================================================
adjons          = TR*nskip;

% | RUNTIME OPTIONS
% | ===========================================================================
if runtest, subdir = subdir(1); end

% | SUBJECT LOOP
% | ===========================================================================
allinput = [];
for s = 1:length(subdir)

    % | Check Subject and Define Folders
    % | ========================================================================
    rundir      = files([subdir{s} filesep 'raw' filesep runid]);
    if isempty(rundir), printmsg('Valid run directory not found, moving on...', 'msgtitle', subnam{s}); continue; end
    analysisdir = fullfile(subdir{s}, 'analysis', analysisname);
    if any([exist(fullfile(analysisdir, 'mask.img'), 'file') exist(fullfile(analysisdir, 'mask.nii'), 'file')])
        printmsg('Level 1 job probably already estimated, moving on...', 'msgtitle', subnam{s}); continue;
    end
    printmsg(sprintf('Building Level 1 Job for %d Runs', length(rundir)),'msgtitle', subnam{s});

    % | Behavioral and Nuisance Regressor Files
    % | ========================================================================
    nuisance    = files([subdir{s} filesep 'raw' filesep runid filesep nuisancepat]);
    behav       = files([subdir{s} filesep 'behav' filesep behavid]);

    % | Get Images
    % | ========================================================================
    images          = cell(size(rundir));
    if ~isempty(fnepi)
        subidx = strcmp(fnepi.subname, subnam{s});
        images = fnepi.epifname(subidx);
    else
        for r = 1:length(rundir)
            images{r} = files([rundir{r} filesep epipat]);
            if isempty(images{r})
                error('\nImage data not found! Failed search pattern:\n%s', [rundir{r} filesep epipat]);
            end
        end
    end

    % | Run Loop
    % | ========================================================================
    for r = 1:length(rundir)

        b = get_behavior(behav{r});
        b.blockwise(:,3) = b.blockwise(:,3) - adjons;

        % | Columns for b.blockwise
        % | =====================================================================
        % 1 - Block
        % 2 - Cond
        % 3 - Onset
        % 4 - Duration
        % 5 - Total_Errors
        % 6 - Foil_Errors

        % | Conditions
        % | =====================================================================
        for c = 1:length(b.condlabels)
            runs(r).conditions(c).name      = b.condlabels{c};
            runs(r).conditions(c).onsets    = b.blockwise(b.blockwise(:,2)==c, 3);
            runs(r).conditions(c).durations = b.blockwise(b.blockwise(:,2)==c, 4);
        end

        % | Floating Parametric Modulators
        % | =====================================================================
        if ~isempty(covidx)
            allpm           = b.blockwise(:,4:6);
            modelpm         = allpm(:,covidx);
            modelpmnames    = pmnames;
            novaridx = find(nanstd(modelpm)==0);
            if ~isempty(novaridx), modelpm(:,novaridx) = []; modelpmnames(novaridx) = []; end
            for p = 1:length(modelpmnames)
                runs(r).floatingpm(p).name = modelpmnames{p};
                runs(r).floatingpm(p).onsets = b.blockwise(:,3);
                runs(r).floatingpm(p).durations = b.blockwise(:,4);
                runs(r).floatingpm(p).values = modelpm(:,p);
            end
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
    general_info.mt_onset           = 8;

    % | Contrasts
    % | ========================================================================
    ncond   = length(b.condlabels);
    w1      = eye(ncond);
    w2      = [1 0 -1; 0 1 -1; 1 -1 0; .5 .5 -1];
    weights = [w1; w2];
    ncon    = size(weights,1);
    for c = 1:ncon
        contrasts(c).type       = 'T';
        contrasts(c).weights    = weights(c,:);
        contrasts(c).name       = bspm_conweights2names(weights(c,:), b.condlabels);
    end
    if fcontrast
        contrasts(ncon+1).type      = 'F';
        contrasts(ncon+1).name      = 'Omnibus';
        contrasts(ncon+1).weights   = eye(ncond);
    end

    % | Make Job
    % | ========================================================================
    allinput{s} = bspm_level1(images, general_info, runs, contrasts);

    % | Cleanup Workspace
    % | ========================================================================
    clear general_info runs contrasts b modelpm modelpmnames images

end
end
% =========================================================================
% * SUBFUNCTIONS
% =========================================================================
function b = get_behavior(in)
% GET_BEHAVIOR
%
%   USAGE: b = get_behavior(in)
%
%       in      behavioral data filename (.mat)
%
%       Columns for b.blockwise
%          1 - Block
%          2 - Cond
%          3 - Onset
%          4 - Duration
%
% CREATED: Bob Spunt, Ph.D. (bobspunt@gmail.com) - 2014.02.24
% =========================================================================
if nargin < 1, error('USAGE: b = get_behavior(in)'); end
if iscell(in), in = char(in); end

% | read data
% | ========================================================================
d = load(in);
b.subjectID = d.subjectID;
if ismember({'result'},fieldnames(d))
    data        = d.result.trialSeeker;
    blockwise   = d.result.blockSeeker;
else
    data        = d.trialSeeker;
    blockwise   = d.blockSeeker;
end
blockwise(:,3:4) = 0;

% | blockwise onsets and durations
% | ========================================================================
nblocks     = size(blockwise, 1);
for i = 1:nblocks
   cblock = data(data(:,1)==i,:);
   blockwise(i, 3) = cblock(1, 6);
   blockwise(i, 4) = cblock(end, 9) - cblock(1, 6);
end
b.blockwise = blockwise;
b.condlabels = {'Face' 'Hand' 'Scram'};
b.varlabels = {'Block' 'Cond' 'Onset' 'Duration'};
b.num_hits = sum(data(:,3)==4 & data(:,7)>0);
b.num_misses = sum(data(:,3)~=4 & data(:,7)>0);
end
function mfile_showhelp(varargin)
% MFILE_SHOWHELP
ST = dbstack('-completenames');
if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
eval(sprintf('help %s', ST(2).file));
end

