function allinput = wrapper_level1_peers_conte_loi2(covidx, varargin)
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
            'HPF',              100,            ...
            'armethod',         2,              ...
            'epifname',         [],             ...
            'nuisancepat',      'bad*txt',      ...
            'epipat',           'swbu*nii',      ...
            'subid',            'RA*',          ...
            'runid',            'EP*LOI_2*',     ...
            'behavid',          'loi2*mat',     ...
            'basename',         'LOI2',         ...
            'brainmask',        '',             ...
            'model',            '2X2',          ...
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
studydir = '/Users/bobspunt/Documents/fmri/peers';
if strfind(pwd,'/home/spunt'), studydir = '/home/spunt/data/peers'; end
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
analysisname  = sprintf('%s_%s_Pmodby_%s_%s_%ds_%s', basename, model, ...
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
    if ncond==2
        w2 = [1 -1]; 
    else
        w2 = [.5 .5 -.5 -.5; .5 -5 .5 -.5; 1 0 -1 0; 0 1 0 -1; 1 -1 0 0; 0 0 1 -1];
    end
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

% | re-code data
% | ========================================================================
con     = blockwise(:,2);
condata = data(:,3); 
switch lower(opt)
    case {'2x2'}
        b.condlabels = {'Why_Face' 'Why_Hand' 'How_Face' 'How_Hand'};
    case {'1x2'}
        b.condlabels = {'Why' 'How'};
        blockwise(con==2, 2)                    = 1;
        data(condata==2, 3)                     = 1; 
        blockwise(ismember(con,[3 4]), 2)       = 2;
        data(ismember(condata,[3 4]), 3)        = 2; 
end
for i = 1:length(unique(data(:,3)))
    cdata = data(data(:,3)==i, [7 10]);
    b.accuracy(i)  = 100*(sum(cdata(:,2)==0)/size(cdata,1));
    b.rt(i)        = nanmean(cdata(:,1)); 
end
b.blockwise = blockwise;
b.varlabels = {'Block' 'Cond' 'Onset' 'Duration' 'Total_Errors' 'Foil_Errors'};  
end    
function mfile_showhelp(varargin)
% MFILE_SHOWHELP
ST = dbstack('-completenames');
if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
eval(sprintf('help %s', ST(2).file));  
end

