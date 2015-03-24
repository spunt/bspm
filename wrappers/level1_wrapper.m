function matlabbatch = level1_wrapper(covidx,varargin)
% matlabbatch = level1_wrapper(covidx,varargin)
%
% To show default settings, run without any arguments.
%
% COVIDX 
%       01 - Duration          
%       02 - Errors (Total)
%       03 - Errors (Foils)
% 

% | SET DEFAULTS AND PARSE VARARGIN
% | ===========================================================================
defaults = {
            'HPF',              128,            ...
            'armethod',         1,              ... 
            'nuisancepat',      '',             ...
            'epipat',           'lswbua*nii',   ...
            'subid',            'RA*',          ...
            'runid',            'EP*LOIS*',     ...
            'behavid',          'loi*mat',      ...
            'basename',         'LOIS',         ...
            'brainmask',        '',             ...
            'model',            '2x3',          ...
            'fcontrast',        0,              ...
            'nskip',            2,              ...
            'runtest',          0,              ...
            'is4D',             0               ...
             };
vals = setargs(defaults, varargin);
if nargin==0, mfile_showhelp; fprintf('\t= DEFAULT SETTINGS =\n'); disp(vals); return; end
fprintf('\n\t= CURRENT SETTINGS =\n'); disp(vals); 

% | PATHS
% | ===========================================================================
studydir = '/Users/bobspunt/Documents/fmri/lois';
if strfind(pwd,'/home/spunt'), studydir = '/home/spunt/data/lois'; end
[subdir, subnam] = files([studydir filesep subid]);

% | ANALYSIS NAME
% | ===========================================================================
armethodlabels  = {'NoAR1' 'AR1' 'WLS'};
covnames        = {'Duration' 'Errors' 'FoilErrors'};
pmnames         = regexprep(covnames(covidx), '_', '');
pmstr           = sprintf(repmat('_%s', 1, length(pmnames)), pmnames{:}); pmstr(1)= [];
analysisname    = sprintf('%s_%s_Pmodby_%s_%s_%ds_%s', basename, model, ...
                        pmstr, armethodlabels{armethod + 1}, HPF, bob_timestamp);
printmsg(analysisname, 'Analysis Name'); 
                                    
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
    if isempty(rundir), printmsg('Valid run directory not found, moving on...', subnam{s}); continue; end
    analysisdir = fullfile(subdir{s}, 'analysis', analysisname); 
    if any([exist(fullfile(analysisdir, 'mask.img'), 'file') exist(fullfile(analysisdir, 'mask.nii'), 'file')])
        printmsg('Level 1 job probably already estimated, moving on...', subnam{s}); continue; 
    end
    printmsg(sprintf('Building Level 1 Job for %d Runs', length(rundir)), subnam{s}); 

    % | Behavioral and Nuisance Regressor Files
    % | ========================================================================
    nuisance    = files([subdir{s} filesep 'raw' filesep runid filesep nuisancepat]);
    behav           = files([subdir{s} filesep 'behav' filesep behavid]);
    
    % | Run-Level Parameters
    % | ========================================================================
    epihdr          = files([rundir{1} filesep 'dicom*mat']);
    epidicominfo    = bspm_get_dicom_info(char(epihdr), 0);
    TR              = epidicominfo.parameterinfo.TR/1000;
    adjons          = TR*nskip;
    mt_res          = epidicominfo.sliceinfo.number;
    mt_onset        = round(mt_res/2);

    % | Run Loop
    % | ========================================================================
    images          = cell(size(rundir)); 
    for r = 1:length(rundir)

        % | Data for Current Run
        % | =====================================================================
        images{r}   = files([rundir{r} filesep epipat]);
        b = get_behavior(behav{r}, model);
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
    general_info.mt_res             = mt_res;
    general_info.mt_onset           = mt_onset;
    general_info.hrf_derivs         = [0 0]; 

    % | Contrasts
    % | ========================================================================
    ncond   = length(b.condlabels);
    w1      = eye(ncond);
    w2      = [ 1 0 0 -1 0 0; 0 1 0 0 -1 0; 0 0 1 0 0 -1; -.5 .5 0 .5 -.5 0; ...
                -.5 0 .5 .5 0 -.5; -.5 .5 0 -.5 .5 0; -.5 0 .5 -.5 0 .5];
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
    matlabbatch = [matlabbatch bspm_level1(images, general_info, runs, contrasts)]; 

    % | Cleanup Workspace
    % | ========================================================================
    clear general_info runs contrasts b modelpm modelpmnames

end
end
% =========================================================================
% * SUBFUNCTIONS
% =========================================================================
function b = get_behavior(in, opt)
% GET_BEHAVIOR
%
%   USAGE: b = get_behavior(in, opt)
%       
%       in      behavioral data filename (.mat)
%       opt     '2x3to2x2'  - social vs. nonsocial
%               '2x3'       - face vs. hand vs. nonsocial
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
if nargin < 1, error('USAGE: b = get_behavior(in, opt)'); end
if nargin < 2, opt = '2x3'; end
if iscell(in), in = char(in); end
 
% | read data
% | ========================================================================
d = load(in);
b.subjectID = d.subjectID; 
if ismember({'result'},fieldnames(d))
    data        = d.result.trialSeeker;
    blockwise   = d.result.blockSeeker; 
    items       = d.result.preblockcues(d.result.blockSeeker(:,4));
else
    data        = d.trialSeeker;
    blockwise   = d.blockSeeker; 
    items       = d.ordered_questions;
end

% | blockwise accuracy and durations
% | ========================================================================
ntrials         = length(data(data(:,1)==1,1));
data(:,10)      = data(:,4)~=data(:,8); % errors
data(data(:,8)==0, 7:8) = NaN; % NR to NaN
data(:,9)       = data(:,6) + data(:,7); 
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
con = blockwise(:,2);
switch opt
    case {'2x3'}
        b.condlabels = {'Why_NS' 'Why_Face' 'Why_Hand' 'How_NS' 'How_Face' 'How_Hand'};
    case {'2x3to2x2'}
        b.condlabels = {'WhyFace' 'WhyHand' 'HowFace' 'HowHand'};
        blockwise(ismember(con,[2 3]), 2)       = 2;
        blockwise(con==4, 2)                    = 3;
        blockwise(ismember(con,[5 6]), 2)       = 4;
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

