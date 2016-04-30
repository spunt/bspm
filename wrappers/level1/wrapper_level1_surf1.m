function allinput = wrapper_level1_surf1(pmidx, varargin)
% mallinput = wrapper_level1_surf1(pmidx, varargin)
%
% To show default settings, run without any arguments.
%
%     PMIDX 
%          1 - Signed Valence Rating (1=Very Bad, 9=Very Good)
%          2 - Signed Valence Rating RT
%          3 - Unsigned Valence Rating 
%          4 - Understanding Rating (1=Not at all, 9=Completely)
%          5 - Understanding Rating RT
%  

% | SET DEFAULTS AND PARSE VARARGIN
% | ===========================================================================
defaults = {
         'studydir',    '/Users/bobspunt/Documents/fmri/dog', ...
         'studyname',   'dog',                                ...
         'HPF',         100,                                  ...
         'armethod',    2,                                    ...
         'nuisancepat', 'rp*txt',                             ...
         'epipat',      'swbua*nii*',                         ...
         'subid',       'RA*',                                ...
         'runid',       'EP*SURF1*',                          ...
         'tag',         's6w2rp',                             ...
         'behavid',     'surf1*mat',                          ...
         'rateid',      'rate*mat',                           ...
         'basename',    'SURF1',                              ...
         'brainmask',   '',                                   ...
         'fcontrast',   0,                                    ...
         'nskip',       4,                                    ...
         'runtest',     0,                                    ...
         'is4D',        1,                                    ...
         'TR',          1,                                    ...
         'yesnokeys',   [1 2]                                 ...
         };
vals = setargs(defaults, varargin);
if nargin==0, mfile_showhelp; fprintf('\t= DEFAULT SETTINGS =\n'); disp(vals); return; end
fprintf('\n\t= CURRENT SETTINGS =\n'); disp(vals); 

% | PATHS
% | ===========================================================================
if strfind(pwd,'/home/spunt'), studydir = fullfile('/home/spunt/data', studyname); end
[subdir, subnam] = files([studydir filesep subid]);

% | ANALYSIS NAME
% | ===========================================================================
armethodlabels  = {'NoAR1' 'AR1' 'WLS'};
covnames        = {'Signed_Valence' 'Valence_RT' 'Unsigned_Valence' 'Understanding' 'Understanding_RT'};
if ~isempty(pmidx)
    pmnames         = regexprep(covnames(pmidx), '_', '');
    pmstr           = sprintf(repmat('_%s', 1, length(pmnames)), pmnames{:}); pmstr(1)= [];
else
    pmstr = 'None'; 
end
analysisname  = sprintf('%s_%s_Pmodby_%s_%s_%ds_%s', basename, tag, ...
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
    printmsg(sprintf('Building Level 1 Job for %d Runs', length(rundir)),'msgtitle', subnam{s}); 

    % | Behavioral and Nuisance Regressor Files
    % | ========================================================================
    nuisance    = files([subdir{s} filesep 'raw' filesep runid filesep nuisancepat]);
    behav       = files([subdir{s} filesep 'behav' filesep behavid]);
    rate        = files([subdir{s} filesep 'behav' filesep rateid]); 
    
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
        b = get_behavior(behav{r}, rate);
        b.data(:,3) = b.data(:,3) - adjons;
        
        % | Columns for b.data
        % | =====================================================================
        %  1 - Trial #
        %  2 - Cond
        %  3 - Onset
        %  4 - Duration
        %  5 - Signed Valence Rating (1=Very Bad, 9=Very Good)
        %  6 - Signed Valence Rating RT
        %  7 - Unsigned Valence Rating 
        %  8 - Understanding Rating (1=Not at all, 9=Completely)
        %  9 - Understanding Rating RT
        
        if ~isempty(pmidx)
            allpm           = b.data(:,5:9);
            modelpm         = allpm(:,pmidx);
            modelpmnames    = pmnames; 
            novaridx = find(nanstd(modelpm)==0);
            if ~isempty(novaridx), modelpm(:,novaridx) = []; modelpmnames(novaridx) = []; end
        end
        
        % | Conditions
        % | =====================================================================
        allcondname = []; 
        for c = 1:length(b.condlabels)
            runs(r).conditions(c).name      = b.condlabels{c};
            allcondname = [allcondname {b.condlabels{c}}]; 
            runs(r).conditions(c).onsets    = b.data(b.data(:,2)==c, 3); 
            runs(r).conditions(c).durations = b.data(b.data(:,2)==c, 4);
            if ~isempty(pmidx) && c < 4
                for p = 1:length(modelpmnames)
                    runs(r).conditions(c).parameters(p).name    = strcat(b.condlabels{c}, '_', modelpmnames{p});
                    allcondname = [allcondname strcat(b.condlabels{c}, {'_'}, modelpmnames{p})]; 
                    runs(r).conditions(c).parameters(p).values  = modelpm(b.data(:,2)==c, p) - mean(modelpm(b.data(:,2)==c, p));
                end
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
    general_info.mt_onset           = 1;

    % | Contrasts
    % | ========================================================================
    ncond   = length(allcondname); 
    
    rmidx   = ismember(allcondname, {'Scramble' 'Catch'});
    w1      = eye(ncond);
    w1(find(rmidx), :) = [];
    if ~isempty(pmidx)
        pmpad = repmat(0, 1, length(modelpmnames));
    else
        pmpad = [];
    end
    w2 = [  1 pmpad 0 pmpad 0 pmpad -1 0; 
            0 pmpad 1 pmpad 0 pmpad -1 0; 
            0 pmpad 0 pmpad 1 pmpad -1 0;
            1 pmpad -1 pmpad 0 pmpad 0 0; 
            1 pmpad 0 pmpad -1 pmpad 0 0;
            0 pmpad 1 pmpad -1 pmpad 0 0    ]; 
    weights = [w1; w2]; 
    ncon    = size(weights,1);
    for c = 1:ncon
        contrasts(c).type       = 'T';
        contrasts(c).weights    = weights(c,:);
        contrasts(c).name       = bspm_conweights2names(weights(c,:), allcondname);
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
function b = get_behavior(in, rate)
% GET_BEHAVIOR
%
%   USAGE: b = get_behavior(in)
%       
%       in      behavioral data filename (.mat)
%
%       Columns for b.data
%          1 - Trial #
%          2 - Cond
%          3 - Onset
%          4 - Duration
%          5 - Signed Valence Rating (1=Very Bad, 9=Very Good)
%          6 - Signed Valence Rating RT
%          7 - Unsigned Valence Rating 
%          8 - Understanding Rating (1=Not at all, 9=Completely)
%          9 - Understanding Rating RT
%
% RATINGS
%   1 - How does the photograph make you feel? (1=Very Bad, 9=Very Good)
%   2 - Do you understand what he or she is feeling? (1=Not at all, 9=Completely)
%
% CREATED: Bob Spunt, Ph.D. (bobspunt@gmail.com) - 2014.02.24
% =========================================================================
if nargin < 2, error('USAGE: b = get_behavior(in, rate)'); end
if iscell(in), in = char(in); end
if iscell(rate), rate = char(rate); end

% | read task data
% | ========================================================================
d = load(in);
b.subjectID     = d.subjectID;
[~,b.stimulus]  = cellfun(@fileparts, d.slideName', 'unif', false);
b.condlabels    = {'Human' 'Monkey' 'Dog' 'Scramble' 'Catch'};
b.varlabels     = {'Trial' 'Cond' 'Onset' 'Duration' 'Signed_Valence' 'Valence_RT' 'Unsigned_Valence' 'Understanding' 'Understanding_RT'};
data            = d.Seeker;
b.percentcaught = 100*(sum(data(:,3)==1 & data(:,7)==1)/sum(data(:,3)==1));
data(data(:,3)==1, 2) = 5; % catch trials 
data(data(:,8)==0, 8) = 1.75; % duration
stimidx = data(data(:,2) < 4, 4);
b.data = data(:,[1 2 6 8]);
b.data(:,5:8) = NaN;

% | read rating data
% | ========================================================================
r = load(rate); 
ratedata        = r.Seeker;
[~, ratestim]   = cellfun(@fileparts, r.slideName', 'unif', false);

% 1 - How does the photograph make you feel? (1=Very Bad, 9=Very Good)
valence         = ratedata(ratedata(:,2)==1, [1 3 4]); 
valence         = sortrows(valence, 1);
b.data(data(:,2)<4, 5:6) = valence(stimidx, 2:3);
b.data(data(:,2)<4, 7) = abs(5 - valence(stimidx,2));

% 2 - Do you understand what he or she is feeling? (1=Not at all, 9=Completely)
understand      = ratedata(ratedata(:,2)==2, [1 3 4]);
understand      = sortrows(understand, 1); 
b.data(data(:,2)<4, 8:9) = understand(stimidx, 2:3); 

end    

