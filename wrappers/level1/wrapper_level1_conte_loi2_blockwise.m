function allinput = wrapper_level1_conte_loi2_blockwise(covidx, varargin)
% matlabbatch = wrapper_level1_conte(covidx, varargin)
%
% To show default settings, run without any arguments.
%
%     COVIDX
%       01 - Duration
%       02 - Errors (Total)
%       03 - Errors (Foils)
%

% | SET DEFAULTS AND PARSE VARARGIN
% | ===========================================================================
defaults = {
            'studydir',         '/Users/bobspunt/Documents/fmri/conte',            ...
            'HPF',              100,            ...
            'armethod',         2,              ...
            'epifname',         [],             ...
            'nuisancepat',      's6w3_badscan*txt',      ...
            'epipat',           's6w3bua*nii',  ...
            'subid',            'CC*',          ...
            'runid',            'EP*LOI_2*',    ...
            'behavid',          'loi2*mat',     ...
            'basename',         'LOI2_s6w3',    ...
            'brainmask',        '',             ...
            'model',            'BLOCKWISE',    ...
            'fcontrast',        0,              ...
            'maskthresh',       0.25,           ...
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
if strfind(pwd,'/home/spunt'), studydir = '/home/spunt/data/conte'; end
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

analysisname  = sprintf('%s_%s_Pmodby_%s_%s_%ds_ImpT%d_%s', basename, model, ...
                        pmstr, armethodlabels{armethod + 1}, HPF, maskthresh*100, bob_timestamp);
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

        % | Data for Current Run
        % | =====================================================================
        b = get_behavior(behav{r}, model, yesnokeys);
        b.blockwise(:,3) = b.blockwise(:,3) - adjons;
        
        % | Sort by condlabel so betas refer to same block for all subjects
        % | =====================================================================
        data = [b.condlabels num2cell(b.blockwise)];
        data = sortrows(data, -1);
        b.condlabels = data(:,1);
        b.blockwise = cell2mat(data(:,2:end));

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
            runs(r).conditions(c).onsets    = b.blockwise(c, 3);
            runs(r).conditions(c).durations = b.blockwise(c, 4);
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
    general_info.maskthresh         = maskthresh;
    general_info.hrf_derivs         = [0 0];
    general_info.mt_res             = 16;
    general_info.mt_onset           = 8;

    % | Contrasts
    % | ========================================================================
    ncond   = length(b.condlabels);
    cond    = b.blockwise(:,2)';
    w1      = (ismember(cond, [1 2]) - ismember(cond, [3 4]))/2;
    w2      = (ismember(cond, [1]) - ismember(cond, [3]))/2;
    w3      = (ismember(cond, [2]) - ismember(cond, [4]))/2;
    weights = [w1; w2; w3];
    ncon    = size(weights,1);
    conname = {'Why_-_How' 'Face_Why_-_How' 'Hand_Why_-_How'};
    for c = 1:ncon
        contrasts(c).type       = 'T';
        contrasts(c).weights    = weights(c,:);
        contrasts(c).name       = conname{c};
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
    clear general_info runs contrasts b modelpm modelpmnames

end
end
% =========================================================================
% * SUBFUNCTIONS
% =========================================================================
function b = get_behavior(in, opt, yesnokeys)
% GET_BEHAVIOR
%
%   USAGE: b = get_behavior(in, opt)
%
%       in      behavioral data filename (.mat)
%       opt     '2x2'  - full design
%               '1x2'  - why vs. how
%
%       Columns for b.blockwise
%          1 - Block
%          2 - Cond
%          3 - Onset
%          4 - Duration
%          5 - Total_Errors
%          6 - Foil_Errors
%
% CREATED: Bob Spunt, Ph.D. (bobspunt@gmail.com) - 2014.02.24
% =========================================================================
if nargin < 1, error('USAGE: b = get_behavior(in, opt, yesnokeys)'); end
if nargin < 2, opt = '2x2'; end
if nargin < 3, yesnokeys = [1 2]; end
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
qidx = blockwise(:, end);
questions = regexprep(d.result.preblockcues(qidx), 'Is the person ', '');
questions = regexprep(questions, ' ', '_');
questions = regexprep(questions, '?', '');
questions = regexprep(questions, '-', '_'); 

% | blockwise accuracy and durations
% | ========================================================================
ntrials         = length(data(data(:,1)==1,1));
data(data(:,8)==yesnokeys(1), 8) = 1;
data(data(:,8)==yesnokeys(2), 8) = 2;
data(:,10)      = data(:,4)~=data(:,8); % errors
data(data(:,8)==0, 7:8) = NaN; % NR to NaN
blockwise(:,3)  = data(data(:,2)==1, 6);
blockwise(:,4)  = data(data(:,2)==ntrials, 9) - data(data(:,2)==1, 6);

% | compute block-wise error counts
% | ========================================================================
for i = 1:size(blockwise, 1)
    blockwise(i,5) = sum(data(data(:,1)==i,10));  % all errors
    blockwise(i,6) = sum(data(data(:,1)==i & data(:,4)==2, 10)); % foil errors
end
for i = 1:length(unique(data(:,3)))
    cdata = data(data(:,3)==i, [7 10]);
    b.accuracy(i)  = 100*(sum(cdata(:,2)==0)/size(cdata,1));
    b.rt(i)        = nanmean(cdata(:,1));
end

% | Blockwise Labels
% | ========================================================================
condlabels = {'Why-Face' 'Why-Hand' 'How-Face' 'How-Hand'};
qcond = condlabels(blockwise(:,2))';
rt = round(blockwise(:,4)*1000);
err = blockwise(:,5);
b.condlabels = strcat(upper(qcond), '-', upper(questions), '-', num2str(rt), 'ms', '-', num2str(err), 'error');
b.condlabels = regexprep(b.condlabels, ' ', '');
b.blockwise = blockwise;
b.varlabels = {'Block' 'Cond' 'Onset' 'Duration' 'Total_Errors' 'Foil_Errors'};
end
function mfile_showhelp(varargin)
% MFILE_SHOWHELP
ST = dbstack('-completenames');
if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
eval(sprintf('help %s', ST(2).file));
end

