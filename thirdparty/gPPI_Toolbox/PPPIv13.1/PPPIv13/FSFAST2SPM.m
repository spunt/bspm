function config=FSFAST2SPM(config)
%% This function will convert your FSFAST Analysis to SPM
% Assumes FSFAST datastructure with 3 directories.
%
% Written by Donald McLaren (mclaren.donald@gmail.com)
% GRECC, Edith Nourse Rogers Memorial Veterans Hospital, Bedford, MA
% Department of Neurology, Massachusetts General Hospital, Boston, MA
% Athinoula A. Martinos Center for Biomedical Imaging, Department of Radiology, Charlestown, MA
% Harvard Medical School, Boston, MA
% 03/20/2014

%Example config settings
%config.Study='/autofs/cluster/iaslab/NMASA/mri/func/encoding';
%config.Subject='nmasa_001_121009';
%config.Subject_Subdirectory='bold';
%config.Model='NMASA_N10_neg_encode_misshit';

%% Function begins here
if exist('config','file')==2
    config=load(config);
end
if isstruct(config)
else
    error('Config is not a file, a structure, or a file containing a structure.')
end
while numel(fields(config))==1
    F=fieldnames(config);
    config=config.(F{1}); %Ignore coding error flag.
end

%% Check Config for Correct Parameters
err=0;
if ~isfield(config,'Study')
    disp('Config structure does not contain Study field - study file path')
    err=err+1;
end
if ~isfield(config,'Subject')
    disp('Config structure does not contain Subject field - subject name')
    err=err+1;
end
if ~isfield(config,'Subject_Subdirectory')
    disp('WARNING: Analysis Folder is in the main subject directory, this is unlikely.')
    config.Subject_Subdirectory=[];
end
if ~isfield(config,'Model')
    disp('Config structure does not contain Model field - Name of the model/analysis directory (w/o lh, rh, mni305)')
    err=err+1;
end

if err>0
    disp('Progam will now exit, config file not correct.')
    return
end

if ~exist([config.Study filesep config.Subject filesep config.Subject_Subdirectory],'dir')
    disp('Subject directory not found')
    return
end
if ~exist([config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.Model '.lh'],'dir')
    disp(['Model directory: ' config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.Model '.lh not found'])
    err=err+1;
end
if ~exist([config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.Model '.rh'],'dir')
    disp(['Model directory: ' config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.Model '.rh not found'])
    err=err+1;
end
if ~exist([config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.Model '.mni305'],'dir')
    disp(['Model directory: ' config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.Model '.mni305 not found'])
    err=err+1;
end

if err>0
    return
end

%% Conversion Program Starts Here
model={'.lh' '.rh' '.mni305'};

for mm=1:3
    %% Load model parameters
    config.params=load([config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.Model model{mm} filesep 'X.mat']);
    if exist([config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.params.runflac(1).flac.runlist(1,:) filesep config.params.runflac(1).flac.funcstem '.nii.gz'],'file')
        config.gzip=1;
    elseif exist([config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.params.runflac(1).flac.runlist(1,:) filesep config.params.runflac(1).flac.funcstem '.nii'],'file')
    else
        error(['Functional File Is Missing: ' config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.params.runflac(1).flac.runlist(1,:) filesep config.params.runflac(1).flac.funcstem '.nii/.nii.gz'])
    end
    
    %% SPM results directory
    SPM.swd=[config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep 'SPMana' filesep config.Model model{mm}];
    
    %% Check if SPMmat has been created.
    if exist(fullfile(SPM.swd,'SPM.mat'),'file') == 2 && exist(fullfile(SPM.swd,'mask.img'),'file') == 2
        config.estimate=0;
        continue;
    else
        config.estimate=1;
    end
    
    %% Check for nscan
    try
        if isnumeric(config.params.runflac(1).flac.nruns)
            SPM.nscan=zeros(1,config.params.runflac(1).flac.nruns);
            for jj=1:config.params.runflac(1).flac.nruns
                if config.gzip==1
                    if ~exist([config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.params.runflac(jj).flac.runlist(jj,:) filesep config.params.runflac(jj).flac.funcstem '.nii'],'file');
                        if ~isfield(config,'zip')
                            config.zip=config.gzip;
                        end
                        try
                            system(sprintf('gunzip -c %s > %s',[config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.params.runflac(jj).flac.runlist(jj,:) filesep config.params.runflac(jj).flac.funcstem '.nii.gz'],[config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.params.runflac(jj).flac.runlist(jj,:) filesep config.params.runflac(jj).flac.funcstem '.nii']));
                        catch
                            if usejava('jvm')
                                try
                                    gunzip([config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.params.runflac(jj).flac.runlist(jj,:) filesep config.params.runflac(jj).flac.funcstem '.nii.gz']);
                                catch
                                    if ~exist([config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.params.runflac(jj).flac.runlist(jj,:) filesep config.params.runflac(jj).flac.funcstem '.nii'],'file')
                                        disp('Please unzip the files manually before processing. gunzip is not a valid command on this system.')
                                        error('Files are not unzipped')
                                    end
                                end
                            else
                                if ~exist([config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.params.runflac(jj).flac.runlist(jj,:) filesep config.params.runflac(jj).flac.funcstem '.nii'],'file')
                                    disp('Please unzip the files manually before processing. gunzip is not a valid command on this system.')
                                    error('Files are not unzipped')
                                end
                            end
                        end
                    else
                        if ~isfield(config,'zip')
                            config.zip=0;
                        end
                    end
                end
                SPM.nscan(jj) = config.params.runflac(jj).flac.ntp;
            end
        else
            invokecatchstatement
        end
    catch
        error('Number of scans must be defined.')
    end
    
    %% Fill in SPM Structure
    for jj = 1:numel(SPM.nscan)
        try
            sessdir=[config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.params.runflac(jj).flac.runlist(jj,:)];
            try
                fid = fopen([sessdir filesep config.params.runflac(jj).flac.parfile],'r');
            catch
                error(['Filename: ' sessdir filesep config.params.runflac(jj).flac.parfile ' does not exist. Please check filename. Did you move your data?']);
            end
            if fid==-1
                error(['Filename: ' sessdir filesep config.params.runflac(jj).flac.parfile ' is not readable. Please check filename.']);
            else
                fileline = 0;
                % read ascii-file line by line
                while fileline >= 0
                    % read a line
                    fileline = fgets(fid);
                    if fileline < 0
                        break
                    elseif isempty(strtrim(fileline))
                        fline(2)=0;
                        fileline=0;
                    else
                        for ii=1:4
                            [fl,fileline]=strtok(fileline);
                            fline(ii)=str2double(fl);
                        end
                        fileline=strtrim(fileline);
                    end
                    if fline(2)>0
                        kk=fline(2); %trial type index in FSFAST starts at 0.
                        try
                            SPM.Sess(jj).U(kk).name{1} = SPM.Sess(jj).U(kk).name{1};
                        catch
                            SPM.Sess(jj).U(kk).name{1}=fileline(~isspace(fileline));
                        end
                        if ~isfield(SPM.Sess(jj).U(kk),'ons'); SPM.Sess(jj).U(kk).ons=[];end
                        if ~isfield(SPM.Sess(jj).U(kk),'dur'); SPM.Sess(jj).U(kk).dur=[];end
                        SPM.Sess(jj).U(kk).ons(end+1) = fline(1);
                        SPM.Sess(jj).U(kk).dur(end+1) = fline(3); %% durations are modeled with a duration
                        %Parametric modulators are implemented as conditions in
                        %FSFAST, not as a modulator in SPM. See
                        %(https://surfer.nmr.mgh.harvard.edu/fswiki/FsFastParametricModulation).
                        %This script changes it to be like SPM (see below)
                        SPM.Sess(jj).U(kk).P.name = SPM.Sess(jj).U(kk).name{1};
                        if ~isfield(SPM.Sess(jj).U(kk),'P') || ~isfield(SPM.Sess(jj).U(kk).P,'P'); SPM.Sess(jj).U(kk).P.P=[];end
                        try
                            SPM.Sess(jj).U(kk).P.P(end+1,1) = fline(4);
                            SPM.Sess(jj).U(kk).P.h = 1;
                        catch %#ok<*CTCH>
                            SPM.Sess(jj).U(kk).P.P(end+1,1) = 1;
                            SPM.Sess(jj).U(kk).P.h = 1;
                        end
                    end
                end
            end
        catch
            error(['Timing parameter file: ' sessdir filesep config.params.runflac(jj).flac.parfile ' is not properly formatted.'])
        end
        for kk=1:numel(SPM.Sess(jj))
            %Remove empty names
            if isempty(SPM.Sess(jj).U(kk).name{1})
                SPM.Sess(jj).U(kk)=[];
            end
        end
        for kk=1:numel(SPM.Sess(jj).U)
            %Parametric Modulators
            parmod=[];
            if all(SPM.Sess(jj).U(kk).P.P==1)
                SPM.Sess(jj).U(kk).P.P=[];
                SPM.Sess(jj).U(kk).P.h=0;
                SPM.Sess(jj).U(kk).P.name='none';
                SPM.Sess(jj).U(kk).P;
            else
                parmod(end+1)=kk;
            end
        end
        for kk=parmod
            for ii=1:numel(SPM.Sess(jj))
                if all(SPM.Sess(jj).U(ii).ons==SPM.Sess(jj).U(kk).ons) && ii~=kk
                    SPM.Sess(jj).U(ii).P.P(:,end+1)=SPM.Sess(jj).U(kk).P.P;
                    SPM.Sess(jj).U(kk).P.P=[];
                    SPM.Sess(jj).U(ii).P.h=SPM.Sess(jj).U(kk).P.h;
                    SPM.Sess(jj).U(kk).P.h=0;
                    if strcmp(SPM.Sess(jj).U(kk).P.name{1},'none')
                        SPM.Sess(jj).U(kk).P.name=[];
                    end
                    SPM.Sess(jj).U(ii).P.name{end+1}=SPM.Sess(jj).U(kk).name{1};
                    SPM.Sess(jj).U(kk).P.name{1}='none';
                end
            end
        end
    end
    
    for jj = 1:numel(SPM.nscan)
        SPM.xX.K(jj).HParam = 1/config.params.runflac(jj).flac.hpfCutoffHz; % No highpass filter
    end
    
    %% Add Nuisance Regressors
        for jj = 1:numel(SPM.nscan)
            try
                SPM.Sess(jj).C.C=config.params.runflac(jj).flac.X(:,config.params.runflac(jj).flac.indnuis);
                SPM.Sess(jj).C.C=SPM.Sess(jj).C.C(:,~all(SPM.Sess(jj).C.C==repmat(mean(SPM.Sess(jj).C.C),size(SPM.Sess(jj).C.C,1),1)));
                for kk=1:size(SPM.Sess(jj).C.C,2)
                    SPM.Sess(jj).C.name{kk} = ['nuis' num2str(kk)];
                end
            catch
                SPM.Sess(jj).C.name = cell(0);
                SPM.Sess(jj).C.C = [];
            end
        end
   
    
    %% Check for TR
    try
        if isnumeric(config.params.runflac(1).flac.TR)
            SPM.xY.RT=config.params.runflac(1).flac.TR;
        else
            invokecatchstatement
        end
    catch
        error('TR is not specified.')
    end
    
    %% Obtaining proper file structure for input into SPM.xY fields
    for jj=1:numel(SPM.nscan)
        [path filename ext]=fileparts(config.params.runflac(jj).flac.mri.fspec);
        if strcmp('.gz',ext)
            config.params.P{jj}=[path filesep filename]; %drops .gz field
        else
            config.params.P{jj}=[path filesep filename ext];
        end
    end
    try
        if numel(config.params.P)==numel(SPM.nscan)
            P=[];
            for ii=1:numel(SPM.nscan)
                for kk=1:SPM.nscan(ii)
                    q = [config.params.P{ii} ',' num2str(kk)];
                    P = strvcat(P,q);
                end
                SPM.Sess(ii).row=sum(SPM.nscan(1:ii-1))+1:1:sum(SPM.nscan(1:ii));
            end
            SPM.xY.P=P;
            if size(SPM.xY.P,1)~=sum(SPM.nscan)
                invokecatch
            end
        else
            invokecatch
        end
    catch
        P=[];
        for ii=1:length(SPM.nscan)
            q = spm_select(SPM.nscan(ii),'image',['Select images for session ' num2str(ii)],{},pwd,'.*',['1:' num2str(SPM.nscan(ii))]);
            P = strvcat(P,q);
            SPM.Sess(ii).row=sum(SPM.nscan(1:ii-1))+1:1:sum(SPM.nscan(1:ii));
        end
        SPM.xY.P=P;
    end
    SPM.xY.VY=spm_vol(SPM.xY.P);
    
    %% Optional Settings
    spm('defaults','FMRI')
    global defaults
    try
        Err=0;
        config.olddefs.stats.fmri.fmri_t=spm_get_defaults('stats.fmri.fmri_t');
        if isnumeric(config.MicrotimeRes)
            defaults.stats.fmri.t=config.MicrotimeRes;
        else
            Err=1;
            invokecatchstatement
        end
    catch
        if ~Err
            defaults.stats.fmri.t=spm_get_defaults('stats.fmri.fmri_t');
        else
            error('MicrotimeRes not specified correctly.')
        end
    end
    try
        config.olddefs.stats.fmri.fmri_t0=spm_get_defaults('stats.fmri.fmri_t0');
        if isnumeric(config.MicrotimeOnset)
            defaults.stats.fmri.t0=config.MicrotimeOnset;
        else
            Err=1;
            invokecatchstatement
        end
    catch
        if ~Err
            defaults.stats.fmri.t0=spm_get_defaults('stats.fmri.fmri_t0');
        else
            error('MicrotimeOnset not specified correctly.')
        end
    end
    %SET UP BF OPTIONS
    SPM.xBF.UNITS = 'secs';
    SPM.xBF.dt = config.params.runflac(1).flac.TR/defaults.stats.fmri.t;
    SPM.xBF.T  = defaults.stats.fmri.t;
    SPM.xBF.T0 = -config.params.runflac(jj).flac.stimulusdelay/SPM.xBF.dt;
    if config.params.runflac(1).flac.ana.gammafit==1
        SPM.xBF.name = 'FSFAST Gamma Function';
        %t=0:params.runflac(1).flac.TR:params.runflac(1).flac.ana.timewindow; %
        %Does not match FSFAST, probably because of the u variable in
        %SPM.Sess.U.
        t=0:SPM.xBF.dt:config.params.runflac(1).flac.ana.timewindow;
        delta=config.params.runflac(1).flac.ana.gamdelay;
        tau=config.params.runflac(1).flac.ana.gamtau;
        alpha=config.params.runflac(1).flac.ana.gamexp;
        if numel(delta)~=1 && numel(tau)~=1 && numel(alpha)~=1; error('Only 1 value of delta. alpha, and tau is allowed.'); end
        h = ( ( ((t - delta)./tau).^alpha) .* exp(-((t - delta)./tau)) );
        h(t<delta) = zeros(1,sum(t<delta));
        SPM.xBF.bf = h/((alpha.^alpha)*exp(-alpha));
        [n,m]=size(SPM.xBF.bf);
        if m>n
            SPM.xBF.bf=SPM.xBF.bf';
        end
        SPM.xBF.dt=config.params.runflac(1).flac.TR/defaults.stats.fmri.t;
        SPM.xBF.length=config.params.runflac(1).flac.ana.timewindow;
        SPM.xBF.order=1;
        
    elseif config.params.runflac(1).flac.ana.spmhrffit==1
        if config.params.runflac(1).flac.ana.nspmhrfderiv==0
            SPM.xBF.name='hrf';
        elseif config.params.runflac(1).flac.ana.nspmhrfderiv==1
            SPM.xBF.name='hrf (with time derivative)';
        elseif config.params.runflac(1).flac.ana.nspmhrfderiv==2
            SPM.xBF.name='hrf (with time and dispersion derivatives)';
        else
            error('derivatives are not set correctly.')
        end
        SPM.xBF = spm_get_bf(SPM.xBF);
    else
        error('firfit is not a valid option at this time.')
    end
    SPM.xBF.Volterra = 1;  %% This is not modeling volterra interactions
    
    %Restore defaults, these were changed by spm_firstlevel_checkstruct
    spm_get_defaults('stats.fmri.fmri_t',config.olddefs.stats.fmri.fmri_t); % Restore old timing
    spm_get_defaults('stats.fmri.fmri_t0',config.olddefs.stats.fmri.fmri_t0); % parameters
    SPM.xGX.iGXcalc = 'None';
    SPM.xVi.form = 'none';
    SPM.xGX.sGXcalc = 'mean voxel value';
    SPM.xGX.sGMsca =  'session specific';
    
    % Change Directory
    try
        cd(SPM.swd)
    catch
        mkdir(SPM.swd)
        cd(SPM.swd)
    end

    % Delete any existing files
    delete beta_00*
    delete ResMS.*
    delete RPV.*
    delete mask.*
   
    % Save SPM
    try
        save([SPM.swd filesep 'SPM.mat'],'SPM')
    catch
        mkdir(SPM.swd)
        save([SPM.swd filesep 'SPM.mat'],'SPM')
    end
    
    disp('estimate_SPM.m')
    SPM = spm_fmri_spm_ui(SPM);
    if mm==1 || mm==2
        SPM.xVol.FWHM=[];
        SPM.xVol.VRpv=[];
        SPM.xVol.R=[];
    elseif mm==3
    end
    SPM.xM.TH(:) = -Inf;  %% disable threshold masking
    SPM=spm_spm(SPM);
    
    %Omnibus Contrast for gPPI
    contrast=defContrasts(SPM,0,-1);
    xCon = spm_FcUtil('Set',contrast.name,contrast.STAT,'c',contrast.c,SPM.xX.xKXs);
    init=length(SPM.xCon);
    if init~=0
        SPM.xCon(init+1) = xCon;
    elseif init==0
        SPM.xCon = xCon;
    else
        triggercatchstatement
    end
    SPM = spm_contrasts(SPM,init+1);
    
    if config.gzip==1 && config.zip==1
        for jj=1:numel(SPM.nscan)
            try
                if exist([config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.params.runflac(jj).flac.runlist(jj,:) filesep config.params.runflac(jj).flac.funcstem '.nii.gz'],'file')
                    delete([config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.params.runflac(jj).flac.runlist(jj,:) filesep config.params.runflac(jj).flac.funcstem '.nii']);
                else
                    system(sprintf('gzip %s',[config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.params.runflac(jj).flac.runlist(jj,:) filesep config.params.runflac(jj).flac.funcstem '.nii']));
                end
            catch
                if usejava('jvm')
                    try
                        gzip([config.Study filesep config.Subject filesep config.Subject_Subdirectory filesep config.params.runflac(jj).flac.runlist(jj,:) filesep config.params.runflac(jj).flac.funcstem '.nii']);
                    catch
                        break
                    end
                else
                    break
                end
            end
        end
    end
    clear SPM
end
