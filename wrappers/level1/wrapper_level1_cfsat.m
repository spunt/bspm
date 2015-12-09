function matlabbatch = wrapper_level1_cfsat(covidx, varargin)
% matlabbatch = wrapper_level1_lois(covidx, varargin)
%
% To show default settings, run without any arguments.
%
%   COVIDX 
%   1 - Duration
%   2 - Total_Errors
%   3 - Foil_Errors
%   4 - NoResponse
%   5 - N_Foils
%   6 - Word_Count
%   7 - Char_Count
% 

% | SET DEFAULTS AND PARSE VARARGIN
% | ===========================================================================
defaults = {
            'studydir',         fullfile(getenv('HOME'),'Bob','Writing','Empirical','ASD','data'), ...
            'HPF',              128,                    ...
            'armethod',         2,                      ... 
            'nuisancepat',      'bad*txt',              ...
            'epipat',           'sw*nii',               ...
            'subid',            'RA*',                  ...
            'runid',            'EP*SAT*',              ...
            'behavid',          'cfsat*mat',            ...
            'basename',         'CFSAT_s8w3',           ...
            'brainmask',        bspm_brainmask,         ...
            'model',            '2x2',                  ...
            'fcontrast',        1,                      ...
            'pmcontrast',       1,                      ...
            'nskip',            2,                      ...
            'TR',               2.5,                    ...
            'runtest',          0,                      ...
            'is4D',             1                       ...
             };
vals = setargs(defaults, varargin);
if nargin==0, mfile_showhelp; fprintf('\t= DEFAULT SETTINGS =\n'); disp(vals); return; end
fprintf('\n\t= CURRENT SETTINGS =\n'); disp(vals); 

% | PATHS
% | ===========================================================================
if strfind(pwd,'/home/spunt'), studydir = '/home/spunt/data/conte'; end
[subdir, subnam] = files([studydir filesep subid]);

% | ANALYSIS NAME
% | ===========================================================================
armethodlabels  = {'NoAR1' 'AR1' 'WLS'};
covnames        =   {             
                    'Duration'    
                    'Total_Errors'
                    'Foil_Errors'
                    'NoResponse'
                    'N_Foils'
                    'Word_Count'
                    'Char_Count' 
                    };  
pmnames         = regexprep(covnames(covidx), '_', '');
pmstr           = sprintf(repmat('_%s', 1, length(pmnames)), pmnames{:}); pmstr(1)= [];
analysisname    = sprintf('%s_%s_Pmodby_%s_%s_%ds_%s', basename, model, ...
                    pmstr, armethodlabels{armethod + 1}, HPF, bob_timestamp);
printmsg(analysisname, 'msgtitle', 'Analysis Name');
                                    
% | RUNTIME OPTIONS
% | ===========================================================================           
if runtest, subdir = subdir(1); end

% | SUBJECT LOOP
% | ===========================================================================
matlabbatch = []; 
for s = 1:length(subdir)
    
    % | Check Subject and Define Folders
    rundir      = files([subdir{s} filesep 'raw' filesep runid]);
    if isempty(rundir), printmsg('Valid run directory not found, moving on...', 'msgtitle', subnam{s}); continue; end
    analysisdir = fullfile(subdir{s}, 'analysis', analysisname); 
    if any([exist(fullfile(analysisdir, 'mask.img'), 'file') exist(fullfile(analysisdir, 'mask.nii'), 'file')])
        printmsg('Level 1 job probably already estimated, moving on...', 'msgtitle', subnam{s}); continue; 
    end
    printmsg(sprintf('Building Level 1 Job for %d Runs', length(rundir)), 'msgtitle', subnam{s}); 

    % | Behavioral and Nuisance Regressor Files
    nuisance    = files([subdir{s} filesep 'raw' filesep runid filesep nuisancepat]);
    behav       = files([subdir{s} filesep 'behav' filesep behavid]);
    
    % | Run Loop
    images = cell(size(rundir)); 
    for r = 1:length(rundir)

        % | Data for Current Run
        images{r}   = files([rundir{r} filesep epipat]);
        if isempty(images{r})
            error('\nImage data not found! Failed search pattern:\n%s', [rundir{r} filesep epipat]); 
        end
        b = get_behavior(behav{r});
        b.blockwise(:,3) = b.blockwise(:,3) - (TR*nskip);
        
        % Columns for b.blockwise
        %    1 - Block
        %    2 - Cond
        %    3 - Onset
        %    4 - Duration
        %    5 - # Errors
        %    6 - # Foil Errors
        %    7 - # No Response Trials
        %    8 - # Foils
        %    9 - Average Question Word Count
        %    10 - Average Question Character Count
        
        % | Conditions
        for c = 1:length(b.condlabels)
            runs(r).conditions(c).name      = b.condlabels{c};
            runs(r).conditions(c).onsets    = b.blockwise(b.blockwise(:,2)==c, 3); 
            runs(r).conditions(c).durations = b.blockwise(b.blockwise(:,2)==c, 4);
        end
        
        % | Floating Parametric Modulators
        allpm           = b.blockwise(:,4:end);
        modelpm         = allpm(:,covidx);
        modelpmnames    = pmnames; 
        novaridx        = find(nanstd(modelpm)==0);
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

    % | Contrasts Weights
    ncond   = length(b.condlabels);
    npm     = length(pmnames); 
    tmp     = eye(ncond+npm);
    w1      = tmp(1:ncond,:);
    % | Conditions
    % 1 - How Video
    % 2 - How Text
    % 3 - Why Video
    % 4 - Why Text
    w2 =    [
                -1  -1  1   1  
                -1   1 -1   1  
                -1   0  1   0   
                 0  -1  0   1   
                -1   1  0   0   
                 0   0 -1   1   
                -1   1  1  -1
            ];
    w2pos = w2;w2pos(w2<0) = 0;
    wscale = sum(w2pos, 2);
    w2 = w2./repmat(wscale, 1, size(w2, 2));
    weights = w1; 
    weights(end+1:end+size(w2,1),:) = 0; 
    weights(size(w1,1)+1:end,1:size(w2,2)) = w2;
    if npm
        b.condlabels = [b.condlabels; pmnames];
        wpm = tmp(ncond+1:end,:);
        wpm(strcmp(pmnames, 'NoResponse'),:) = []; 
        weights = [weights; wpm]; 
    end
    ncon    = size(weights,1);
    % | T-Contrast
    for c = 1:ncon
        contrasts(c).type       = 'T';
        contrasts(c).weights    = weights(c,:);
        contrasts(c).name       = bspm_conweights2names(weights(c,:), b.condlabels);
    end
    
    % | F-Contrast
    if fcontrast
        contrasts(ncon+1).type      = 'F';
        contrasts(ncon+1).name      = 'Omnibus';
        contrasts(ncon+1).weights   = eye(length(b.condlabels));
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
%          5 - # Errors
%          6 - # Foil Errors
%          7 - # No Response Trials
%          8 - # Foils
%          9 - Average Question Word Count
%          10 - Average Question Character Count
%
% CREATED: Bob Spunt, Ph.D. (bobspunt@gmail.com) - 2014.02.24
% =========================================================================
if nargin < 1, error('USAGE: b = get_behavior(in, opt)'); end
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
data = sortrows(data, [1 2]); % sort data by block #, then trial #

% | Conditions
% 1 - How Video
% 2 - How Text
% 3 - Why Video
% 4 - Why Text

% | blockSeeker
%  1 - block #
%  2 - condition 
%  3 - intended block onset
%  4 - actual block onset
%  5 - actual block duration
%  6 - # correct

% | trialSeeker
%  1 - block #
%  2 - trial #
%  3 - cond # 
%  4 - correct/normative answer (0 = No, 1 = Yes)
%  5 - character count
%  6 - word count
%  7 - actual answer
%  8 - response time (s)
%  9 - actual trial onset (question)
%  10 - actual trial offset (action stimulus)


% | blockwise accuracy and durations
% | ========================================================================
ntrials                 = length(data(data(:,1)==1,1)); % ntrials/block
data(data(:,4)==0,4)    = 2;                    % recode: 1=corr, 2=incorr
data(:,11)              = data(:,4)~=data(:,7); % errors
data(data(:,7)==0, 7:8) = NaN;                  % NR to NaN
blockwise(:,3)          = data(data(:,2)==1,9); % block onsets
blockwise(:,4)          = data(data(:,2)==ntrials,10) - data(data(:,2)==1,9); % block durations

% | compute block-wise error counts
% | ========================================================================
for i = 1:size(blockwise, 1)
    blockwise(i,5) = sum(data(data(:,1)==i,11));  % all errors
    blockwise(i,6) = sum(data(data(:,1)==i & data(:,4)==2, 11)); % foil errors
    blockwise(i,7) = sum(isnan(data(data(:,1)==i,7))); % no response
    blockwise(i,8) = sum(data(:,1)==i & data(:,4)==2); % n foils
    blockwise(i,9) = mean(data(data(:,1)==i,4)); % mean word count
    blockwise(i,10) = mean(data(data(:,1)==i,5)); % mean char count
end

% | re-code data
% | ========================================================================
b.blockwise  = blockwise;
b.condlabels = {
              'How_Video' 
              'How_Text'
              'Why_Video'
              'Why_Text' 
              };
b.varlabels = {              
             'Block'       
             'Cond'        
             'Onset'       
             'Duration'    
             'Total_Errors'
             'Foil_Errors'
             'NoResponse'
             'N_Foils'
             'Word_Count'
             'Char_Count' 
             };           
end    
function mfile_showhelp(varargin)
% MFILE_SHOWHELP
ST = dbstack('-completenames');
if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
eval(sprintf('help %s', ST(2).file));  
end

 
 
 
 
 
 
 
 
