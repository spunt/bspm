function allinput = wrapper_level1_surf2_qwise(covidx, varargin)
% allinput = wrapper_level1_surf2(covidx, varargin)
%
% To show default settings, run without any arguments
%
%     COVIDX
%       01 - RT
%       02 - Errors
%       03 - Foil
%

% | SET DEFAULTS AND PARSE VARARGIN
% | ===========================================================================
defaults = {
           'studydir',    '/Users/bobspunt/Documents/fmri/dog', ...
           'studyname',   'dog',                                ...
           'nuisancepat', 'rp*.txt',                ...
           'epipat',      'swbua*nii',                           ...
           'subid',       'RA*',                                ...
           'runid',       'EP*SURF2*',                          ...
           'behavid',     'surf2*mat',                          ...
           'rateid',      'rate*mat',                           ...
           'basename',    'SURF2_QWISE',                              ...
           'tag',         's6w2rp',                               ...
           'brainmask',   '',                                   ...
           'epifname',    [],                                   ...
           'model',       'QWISE',                                ...
           'HPF',         100,                                  ...
           'armethod',    2,                                    ...
           'junkerrors',  0,                                    ...
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
covnames        = {'RT' 'Err' 'Foil'};
labs            = {'incl' 'excl'};
if ~isempty(covidx)
    pmnames         = regexprep(covnames(covidx), '_', '');
    pmstr           = sprintf(repmat('_%s', 1, length(pmnames)), pmnames{:}); pmstr(1)= [];
else
    pmstr = 'None';
end
analysisname  = sprintf('%s_%s_Pmodby_%s_%sERR_%s_%ds_%s', basename, tag, ...
                        pmstr, labs{junkerrors + 1}, armethodlabels{armethod + 1}, HPF, bob_timestamp);
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
    rate        = files([subdir{s} filesep 'behav' filesep rateid]); 

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
        images{r}   = files([rundir{r} filesep epipat]);
        if isempty(images{r}), error('\nImage data not found! Failed search pattern:\n%s', [rundir{r} filesep epipat]); end
        b = get_behavior(behav{r}, rate, yesnokeys, junkerrors);
        b.data(:,4) = b.data(:,4) - adjons;
        ncond = length(b.qwise.condlabels);
        if any(b.data(:,2)==7)
            b.qwise.condidx(b.data(:,2)==7) = ncond + 1;
            b.qwise.condlabels{end + 1} = 'No_Response';
            ncond = ncond + 1; 
        end
            
        % | Columns for b.data
        % | =====================================================================
        % 01 - Trial
        % 02 - Cond
        % 03 - Foil
        % 04 - Onset
        % 05 - Duration
        % 06 - Error (0=No, 1=Yes)
        % 07 - Signed_Valence (1=Very Bad, 9=Very Good)
        % 08 - Valence_RT
        % 09 - Unsigned_Valence 
        % 10 - Understanding (1=Not at all, 9=Completely)
        % 11 - Understanding_RT

        % | Conditions
        % | =====================================================================
        for c = 1:ncond
            runs(r).conditions(c).name      = b.qwise.condlabels{c};
            runs(r).conditions(c).onsets    = b.data(b.qwise.condidx==c, 4);
            runs(r).conditions(c).durations = b.data(b.qwise.condidx==c, 5);
        end

        % | Floating Parametric Modulators
        % | =====================================================================
        if ~isempty(covidx)
            valididx        = b.data(:,2) < 7;
            allpm           = b.data(valididx, [5 6 3]);
            modelpm         = allpm(:,covidx);
            modelpmnames    = pmnames;
            novaridx = find(nanstd(modelpm)==0);
            if ~isempty(novaridx), modelpm(:,novaridx) = []; modelpmnames(novaridx) = []; end
            for p = 1:length(modelpmnames)
                runs(r).floatingpm(p).name = modelpmnames{p};
                runs(r).floatingpm(p).onsets = b.data(valididx, 4);
                runs(r).floatingpm(p).durations = b.data(valididx, 5);
                runs(r).floatingpm(p).values = modelpm(:,p);
            end % for p = 1:length(modelpmnames)
        end % if ~isempty(covidx)

    end % for r = 1:length(rundir)
    
    % | Cleanup for Single Run Models
    % | ========================================================================
    if length(rundir)==1
        images = images{1};
        if iscell(nuisance)
            nuisance = nuisance{1};
        end
    end

    % | General Model Information
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
    % 1 - Why_Human
    % 2 - Why_Primate
    % 3 - Why_Dog
    % 4 - How_Human
    % 5 - How_Primate
    % 6 - How_Dog
    clab = b.qwise.condlabels; 
    ncond = length(clab);
    w0  = zeros(1, ncond);
    why = cellstrfind(clab, 'WHY');
    how = cellstrfind(clab, 'HOW');
    weights = repmat(w0, length(why), 1);
    weights(:,how) = -1/length(how);
    weights(:,why) = eye(length(why));
    conname = regexprep(clab(why), 'WHY-', '');
    conname = strcat(conname, '_-_', 'HOW');
    ncon    = size(weights,1);
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
function b = get_behavior(in, rate, yesnokeys, junkerrors)
    % GET_BEHAVIOR
    % 
    %   USAGE: b = get_behavior(in, rate, yesnokeys, junkerrors)
    % % 
    %       Columns for b.data
    %         01 - Trial
    %         02 - Cond
    %         03 - Foil
    %         04 - Onset
    %         05 - Duration
    %         06 - Error (0=No, 1=Yes)
    %         07 - Signed_Valence (1=Very Bad, 9=Very Good)
    %         08 - Valence_RT
    %         09 - Unsigned_Valence 
    %         10 - Understanding (1=Not at all, 9=Completely)
    %         11 - Understanding_RT
    % 
    %     RATINGS
    %       1 - How does the photograph make you feel? (1=Very Bad, 9=Very Good)
    %       2 - Do you understand what he or she is feeling? (1=Not at all, 9=Completely)
    % 
    % CREATED: Bob Spunt, Ph.D. (bobspunt@gmail.com) - 2014.02.24
    % =========================================================================
    if nargin < 2, error('USAGE: b = get_behavior(in, rate, yesnokeys, junkerrors)'); end
    if nargin < 3, yesnokeys = [1 2]; end
    if nargin < 4, junkerrors = 0; end
    if iscell(in), in = char(in); end
    if iscell(rate), rate = char(rate); end

    % | read task data
    % | ========================================================================
    %% SEEKER column key %%
    % 1 - trial #
    % 2 - condition (1=HH, 2=HM, 3=HD, 4=LH, 5=LM, 6=LD)
    % 3 - correct (normative) response (1=Yes, 2=No)
    % 4 - slide # (corresponds to order in stimulus dir)
    % 5 - question # (corresponds to order in 'qstim', defined in design.mat)
    % 6 - scheduled question onset
    % 7 - scheduled photo onset
    % 8 - actual stimulus onset (s)
    % 9 - actual response [0 if NR]
    % 10 -response time (s) [0 if NR]
    % 11 -(added) normativity

    d = load(in);
    b.subjectID     = d.subjectID;
    [~,b.stimulus]  = cellfun(@fileparts, d.slideName', 'unif', false); 
    b.condlabels    = {'Why_Human' 'Why_Primate' 'Why_Dog' 'How_Human' 'How_Primate' 'How_Dog'};
    b.varlabels     = {'Trial' 'Cond' 'Foil' 'Onset' 'Duration' 'Error' 'Signed_Valence' 'Valence_RT' 'Unsigned_Valence' 'Understanding' 'Understanding_RT'};
    
    % | task performance
    % | ========================================================================
    data                             = d.Seeker;
    ntrial                           = size(data, 1);
    ncond                            = length(unique(data(:,2)));
    data(data(:,9)==yesnokeys(1), 9) = 1;
    data(data(:,9)==yesnokeys(2), 9) = 2;
    data(:,11)                       = data(:,3)~=data(:,9); % errors
    nridx                            = data(:,9)==0; 
    data(nridx, 9:10)                = NaN; % NR to NaN
    data(nridx, 2)                   = 7; 
    data(nridx, 10)                  = 1.5;
    if junkerrors, data(data(:,11)==1, 2) = 7; end
    if any(data(:,2)==7), b.condlabels{7}   = 'Junk'; end
    b.data                           = data(:,[1 2 3 8 10 11]);
    b.data(:,3)                      = b.data(:,3) - 1;
    b.data(:,7:11)                   = NaN;
    
    % | Question Wise
    % | ========================================================================
   stimidx                          = data(:,4);
    qidx                             = data(:,5);
    b.question = d.qstim(data(:,5));
    [~,tmp] = cellfileparts(d.slideName);
    tmp = cellfun(@(x) x(1:2), tmp, 'Unif', false);
    tmp = regexprep(tmp, 'F_', 'Human');
    tmp = regexprep(tmp, 'M_', 'Primate');
    tmp = regexprep(tmp, 'D_', 'Dog');
    qtmp = regexprep(b.question, '\?', '');
    qtmp = caseconvert(qtmp, 'title');
    qtmp = regexprep(qtmp, ' ', '_');
    qtmp = regexprep(qtmp, 'Looking_At_The_Camera', 'Looking_at_Camera'); 
    ctmp = repmat({'WHY'}, size(qtmp));
    ctmp(b.data(:,2)>3) = {'HOW'};
    b.qwise.triallabels = strcat(ctmp, '-', qtmp, '_', tmp');
    [b.qwise.condlabels, ia, b.qwise.condidx] = unique(b.qwise.triallabels); 
    
    % | read rating data
    % | ========================================================================
    
    r = load(rate); 
    ratedata        = r.Seeker;
    [~, ratestim, ext]   = cellfun(@fileparts, r.slideName', 'unif', false);
    b.photo = strcat(ratestim, ext);
    
    % 1 - How does the photograph make you feel? (1=Very Bad, 9=Very Good)
    valence         = ratedata(ratedata(:,2)==1, [1 3 4]); 
    valence         = sortrows(valence, 1);
    b.data(:, 7:8)  = valence(stimidx, 2:3);
    b.data(:, 9)    = abs(5 - valence(stimidx,2));
    
    % 2 - Do you understand what he or she is feeling? (1=Not at all, 9=Completely)
    understand          = ratedata(ratedata(:,2)==2, [1 3 4]);
    understand          = sortrows(understand, 1); 
    b.data(:, 10:11)    = understand(stimidx, 2:3); 
    
    % | summarize
    % | ========================================================================
    b.summary.ntrial        = ntrial;
    b.summary.ncond         = ncond;
    b.summary.percentnoresp = 100*(sum(nridx)/ntrial);
    b.summary.percenterror  = 100*(sum(data(:,11))/ntrial);
    varidx                  = [6 5 7:10];
    b.summary.varlabels     = b.varlabels(varidx);
    b.summary.data          = zeros(ncond, length(varidx));
    for i = 1:ncond
        cdata = b.data(b.data(:,2)==i, varidx);
        tmp = nanmean(cdata);
        tmp(1) = 100*(1-tmp(1));
        b.summary.data(i,:) = tmp; 
    end
      
end    

