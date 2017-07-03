function Params = memMapData(Params)
% This function maps data into disk as blocks and
% creates the corresponding memory pointer objects.
%
% inputs:
% Params - analysis parameters from initParams.m
%
% See also: RUNANALYSIS, INITPARAMS, MEMMAPFILE

% Last updated: 9.11.2010 by Jukka-Pekka Kauppi
% Tampere University of Technology
% Department of Signal Processing
% e-mail: jukka-pekka.kauppi@tut.fi

showTime(1);

Pub = Params.PublicParams;
Priv = Params.PrivateParams;


% check/create paths:
if exist(Pub.dataDestination,'file') ~= 7
    mkdir(Pub.dataDestination)
end
if exist(Priv.subjectDestination,'file') ~= 7
    mkdir(Priv.subjectDestination)
end
if exist(Priv.subjectFiltDestination,'file') ~= 7
    mkdir(Priv.subjectFiltDestination)
end
if exist(Priv.resultsDestination,'file') ~= 7
    mkdir(Priv.resultsDestination)
end
if exist(Priv.statsDestination,'file') ~= 7
    mkdir(Priv.statsDestination)
end
if exist(Priv.PFDestination,'file') ~= 7
    mkdir(Priv.PFDestination)
end
if exist(Priv.withinDestination,'file') ~= 7
    mkdir(Priv.withinDestination)
end
if exist(Priv.phaseDifDestination,'file') ~= 7
    mkdir(Priv.phaseDifDestination)
end
if exist(Priv.PFsessionDestination) ~= 7
    mkdir(Priv.PFsessionDestination)
end

try
    % load pointers:
    %    ParamsNew = Params;
    load([Pub.dataDestination 'memMaps'])
    load([Pub.dataDestination 'Tag'])
    load([Pub.dataDestination Tag])
    disp(['The project with name "' Pub.dataDescription '" already exists in the given directory.'])
    %    ParamsNew = checkUpdates(Params,ParamsNew); % check parameter changes
    %    Pub = ParamsNew.PublicParams;
    %    Priv = ParamsNew.PrivateParams;
    %    if Priv.filterUpdate
    %        disp('No parameter changes specified -> skip memomy mapping')
    %        return
    %    end
    disp(['Updating existing project ' Pub.dataDescription '...'])
catch
    %    disp(lasterr)
    % if pointers do not exist, initialize pointer-struct:
    memMaps.(Priv.origMapName) = [];
    memMaps.(Priv.synchMapName) = [];
    memMaps.(Priv.phaseSynchMapName) = [];
    memMaps.(Priv.filtMapName) = [];
    memMaps.(Priv.withinMapName) = [];
    memMaps.(Priv.cormatMapName) = [];
    memMaps.(Priv.statMapName) = [];
    memMaps.(Priv.phaseMapName) = [];
    memMaps.(Priv.PFMapName) = [];
    memMaps.(Priv.PFmatMapName) = [];
    memMaps.(Priv.PFMapSessionName) = [];
    memMaps.(Priv.PFmatMapSessionName) = [];
 save([Pub.dataDestination 'memMaps'],'memMaps')
    disp(['Establishing a new project with name "' Pub.dataDescription '"...'])
end


if Pub.calcStandard
    % standard analysis
    if Priv.filterUpdate
        % map preprocessed data:
        updt=mapData('subject',Pub,Priv);
        if(updt)
            load([Pub.dataDestination Pub.dataDescription]); %to ensure that possible updates are here
            Pub = Params.PublicParams;
            Priv = Params.PrivateParams;
        end
        % map subband data:
        if Pub.nrFreqBands > 0
            [~]=mapData('subjectFilt',Pub,Priv);
            
            if Pub.freqCompOn
                % map frequency comparison data:            
                % [~]=mapData('PF',Pub,Priv);
                [~]=mapData('PFmats',Pub,Priv);
            end
        end
    end
    % map basic ISC data:
    [~]=mapData('maps',Pub,Priv);
    % map time-varying curves:
    if Pub.useTemplate
        [~]=mapData('syncCurves',Pub,Priv);
    end
end

if Pub.sessionCompOn
    [~]=mapData('PFmatsSession',Pub,Priv);
end
       
if Pub.calcCorMatrices
    % map correlation matrices:
    [~]=mapData('cormats',Pub,Priv);
end

if Pub.calcStats
    % median, quartile, t and variance ISC maps
    [~]=mapData('stats',Pub,Priv);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% recent stuff:
% map within-subject-data:
%if Pub.useTemplate
%    mapData('within',Pub,Priv);
%end

% inter-subject phase synchronization data:
if Pub.calcPhase
    [~]=mapData('phase',Pub,Priv);
    if Pub.useTemplate
        [~]=mapData('phaseSyncCurves',Pub,Priv);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Params.PublicParams = Pub;
Params.PrivateParams = Priv;

% Return information about computer on which MATLAB is running:
[~,maxsize,endian] = computer;
Params.PrivateParams.computerInfo.endian = endian;
Params.PrivateParams.computerInfo.platform = 'str';
Params.PrivateParams.computerInfo.maxsize = maxsize;

% Params.PrivateParams.filterUpdate
save([Pub.dataDestination Pub.dataDescription],'Params')

showTime(0);

function updt = mapData(typeOfData,Pub,Priv)
% This function memory maps the necessary data.
%
% Inputs:
% typeOfData - either 'subject','subjectfilt', 'maps' or
% 'synch' depending which data is mapped. 'subject' refers to
% original preprocessed data, 'subjectfilt' to wavelet filtered
% data, 'maps' to similarity maps, and 'synch' to synchronization
% curves.
%
% Pub and Priv - parameter structures obtained through
% initParams.m.

updt=false;

switch typeOfData
    case 'subject'
        disp('Mapping pre-processed data:')
        for m = 1:Priv.nrSessions
            for k = 1:Priv.nrSubjects
                disp(['Session ' num2str(m) ', Subject ' num2str(k) ':'])
                fullPath = [Priv.subjectDestination ...
                    Priv.prefixSubject num2str(k) ...
                    Priv.prefixSession num2str(m) ...
                    '.bin'];
                if exist(fullPath,'file') == 2
                    if m == Priv.nrSessions && k == Priv.nrSubjects
                        disp('Subject data already initialized...');return
                    else
                        continue
                    end
                end
                flag = 1;
                while flag ~= 0
                    try
                        if flag == 20
                            error('Problems with writing data, memory mapping quitted!!')
                        end
                        
                        if strcmp(Pub.fileFormatSubj,'nii')
                            I = load_nii(Pub.subjectSource{m,k});
                            I = single(I.img);
                            % Check if original header information and
                            % nifti data size are not the same. If not,
                            % update true data size to params-struct.
                            if k == 1
                                if ~isequal(size(I),Pub.dataSize(m,:))
                                    Pub.dataSize(m,:) = size(I);
                                    Priv.dataSize(m,:) = size(I);
                                    Priv.dataSizeMismatch = true;
                                    Params.PublicParams = Pub;
                                    Params.PrivateParams = Priv;
%                                     Params.PublicParams.dataSize(m,:) = size(I);
%                                     Params.PrivateParams.dataSize(m,:) = size(I);
%                                     Params.PrivateParams.dataSizeMismatch = 1;
%                                     Params.PublicParams = Pub;
                                    save([Pub.dataDestination Pub.dataDescription],'Params')
                                    disp('Original header information does not match true data size,')
                                    disp(['parameter field ''dataSize'' updated and parameters saved to directory ' Pub.dataDestination '.'])                                                            
                                    disp(' ')
                                    
                                    updt=true;
                                end
                            end
                        elseif strcmp(Pub.fileFormatSubj,'mat')
                            I = load(Pub.subjectSource{m,k});
                            fiel = fields(I);
                            I = I.(fiel{1});
                            I = single(I);
                        else
                            error('Extension must be ''nii'' or ''mat''')
                        end
                        
                        if ~exist(fullPath,'file')
                            I = permute(I,[4 2 3 1]);
                            fid = fopen(fullPath, 'w');
                            % write data to a binary file:
                            fwrite(fid, I, 'single');clear I;
                            fclose(fid);clear fid;
                        end
                        % perform dynamic memory mapping:
                        memMap.([Priv.prefixSession num2str(m)]).([...
                            Priv.prefixSubject num2str(k)]) = memmapfile(...
                            fullPath,'format',{'single',...
                            [Priv.dataSize(m,4) Priv.dataSize(m,2:3)],...
                            'tyz'});
                        % range check of the values:
                        M = max(memMap.([Priv.prefixSession ...
                            num2str(m)]).([Priv.prefixSubject ...
                            num2str(k)]).Data(round(Priv.dataSize(m,1)/2)).tyz(:,...
                            round(Priv.dataSize(2)/2),round(Priv.dataSize(m,3)/2)));
                        disp(['Max value (x=' num2str(round(Priv.dataSize(m,1)/2)) ',y='...
                            num2str(round(Priv.dataSize(m,2)/2)) ',z=' num2str(round(Priv.dataSize(m,3)/2)) '): ' num2str(M)])
                        if M > 10e20
                            error('Unexpected value found! Possible cause: processing distributed over different (32bit/64bit) memory architectures.')
                        end
                        flag = 0;
                    catch err
                        if flag == 20
                            error('Problems with writing data, memory mapping quitted!!')
                        end
                        disp(err.message)
                        flag = flag + 1;
                    end
                end
            end
        end
        load([Pub.dataDestination 'memMaps'])
        memMaps = setfield(memMaps,Priv.origMapName,memMap);
        save([Pub.dataDestination 'memMaps'],'memMaps')
        disp(' ')
        
    case 'subjectFilt'
        disp('Mapping sub-band data:')
        for s = 1:Priv.nrSessions
            flg = 1;
            for k = 1:Priv.nrSubjects
                for m = 1:Priv.maxScale + 1
                    fullPath = [Priv.subjectFiltDestination ...
                        Priv.prefixSubjectFilt ...
                        num2str(k) '_' Priv.prefixFreqBand ...
                        num2str(m) '_' Priv.prefixSession ...
                        num2str(s) '_' Priv.transformType ...
                        '.bin'];
                    if exist(fullPath,'file') == 2
                        if k == Priv.nrSubjects && m == Priv.maxScale + 1
                            disp('Filtered subject data already initialized...');return
                        else
                            continue
                        end
                    end
                    if flg == 1
                        I = zeros(Priv.dataSize(s,:),'single') ;
                        I = permute(I,[4 2 3 1]);
                        flg = 0;
                    end
                    
                    flag = 1;
                    while flag ~= 0
                        try
                            if flag == 20
                                error('Problems when mapping data, memory mapping quitted!')
                            end
                            if exist(fullPath,'file') == 0
                                fid = fopen(fullPath, 'w');
                                % write data to a binary file:
                                fwrite(fid, I, 'single');
                                fclose(fid);clear fid;
                            end
                            % perform dynamic memory mapping:
                            memMap.([Priv.prefixSession ...
                                num2str(s)]).([Priv.prefixSubjectFilt ...
                                num2str(k)]).([Priv.prefixFreqBand ...
                                num2str(m)]) = memmapfile(fullPath, ...
                                'format',{'single',[Priv.dataSize(...
                                s,4) Priv.dataSize(s,2:3)]...
                                ,'tyz'},'Writable',true);
                            disp(['Session ' num2str(s) ', Subject '...
                                num2str(k) ', Frequency band ' num2str(m)])
                            flag = 0;
                        catch err
                            if flag == 20
                                error('Problems when mapping data, memory mapping quitted!!')
                            end
                            disp(err.message)
                            flag = flag + 1;
                        end
                    end
                end
            end
        end
        clear I;
        load([Pub.dataDestination 'memMaps'])
        memMaps = setfield(memMaps,Priv.filtMapName,memMap);
        save([Pub.dataDestination 'memMaps'],'memMaps')
        disp(' ')
        
    case 'maps'
        flagMap = [false false];
        I = zeros(Priv.dataSize(1,1:3),'single' );
        disp('Mapping inter-subject synchronization maps:')
        for k = 1:length(Priv.simM)
            if ( ( Pub.ssiOn && strcmp(Priv.simM{k},'ssi') ) || ...
                    ( Pub.nmiOn && strcmp(Priv.simM{k},'nmi') ) || ...
                    ( Pub.corOn && strcmp(Priv.simM{k},'cor') ) || ...
                    ( Pub.kenOn && strcmp(Priv.simM{k},'ken') ) )
                for s = 1:Priv.nrSessions
                    Iw = zeros([Priv.dataSize(s,1:3)...
                        Priv.nrTimeIntervals(s)],'single' );
                    for m = 0:Priv.maxScale + 1
                        disp(['Mapping sim.measure ' num2str(k) ...
                            ', session ' num2str(s) ', band ' num2str(m)])
                        
                        fullPath = [Priv.resultsDestination ...
                            Priv.prefixResults '_' ...
                            Priv.simM{k} '_' Priv.prefixFreqBand ...
                            num2str(m) '_' Priv.prefixSession ...
                            num2str(s) '_' Priv.transformType];
                        if Pub.winOn
                            fullPath_win = [fullPath '_win' '.bin'];
                        end
                        fullPath = [fullPath '.bin'];
                        if Pub.winOn
                            if exist(fullPath,'file') == 2 && exist(fullPath_win,'file') == 2
                                if s == Priv.nrSessions && m == (Priv.maxScale+1)
                                    flagMap = true;
                                else
                                    continue
                                end
                            end
                        else
                            if exist(fullPath,'file') == 2
                                if s == Priv.nrSessions && m == (Priv.maxScale+1)
                                    flagMap = true;
                                else
                                    continue
                                end
                            end
                        end
                        if flagMap
                            disp('ISC maps already mapped...');
                            return
                        end
                        flag = 1;
                        while flag ~= 0
                            try
                                if flag == 20
                                   error('Unexpected problem when initializing data!!')                                 
                                end
                                % map whole session data:
                                if ~exist(fullPath,'file')
                                    fid = fopen(fullPath, 'w');
                                    fwrite(fid, I, 'single');
                                    fclose(fid);clear fid;
                                end
                                memMap.('whole').([Priv.prefixFreqBand ...
                                    num2str(m)]).([Priv.prefixSession ...
                                    num2str(s)]).(Priv.simM{k}) = ...
                                    memmapfile(fullPath,'format',...
                                    {'single',Priv.dataSize(s,1:3)...
                                    ,'xyz'},'Writable',true);
                                if Pub.winOn
                                    % map windowed data:
                                    if ~exist(fullPath_win,'file')
                                        fid = fopen(fullPath_win, 'w');
                                        fwrite(fid, Iw, 'single');
                                        fclose(fid);clear fid;
                                    end
                                    if Priv.nrTimeIntervals(s) == 0
                                        memMap.('win') = [];
                                    else
                                        memMap.('win').([Priv.prefixFreqBand ...
                                            num2str(m)]).([Priv.prefixSession ...
                                            num2str(s)]).(Priv.simM{k}) = ...
                                            memmapfile(fullPath_win,'format',...
                                            {'single',Priv.dataSize(s,1:3)...
                                            ,'xyz'},'Writable',true);
                                    end
                                end
                                flag = 0;
                            catch err
                                if flag == 20
                                   error('Unexpected problem when initializing data!!')                                 
                                end
                                disp(err.message)
                                flag = flag + 1;
                            end
                        end
                    end
                end
            end
        end
        clear I;clear Iw;
        load([Pub.dataDestination 'memMaps'])
        memMaps = setfield(memMaps,Priv.resultMapName,memMap);
        save([Pub.dataDestination 'memMaps'],'memMaps')
        disp(' ')
        
        
    case 'within'
        
        I = zeros([Priv.dataSize(1,1:3) Priv.nrSubjects length(Priv.brainRegions)-1 ],'single') ;
        disp('Mapping within-subject maps:')
        for s = 1:Priv.nrSessions
            %              Iw = single( zeros([Priv.dataSize(s,1:3)...
            %              Priv.nrTimeIntervals(s)]) );
            for m = 0:Priv.maxScale + 1
                disp(['Mapping session ' num2str(s) ', band ' num2str(m)])
                
                fullPath = [Priv.withinDestination ...
                    Priv.prefixWithin '_' ...
                    Priv.simM{3} '_' Priv.prefixFreqBand ...
                    num2str(m) '_' Priv.prefixSession ...
                    num2str(s) '_' Priv.transformType];
                
                %                fullPath_win = [fullPath '_win' '.bin'];
                fullPath = [fullPath '.bin'];
                if exist(fullPath,'file') == 2
                    disp('Data already exist, quit memory mapping...');return
                end
                flag = 1;
                while flag ~= 0
                    try
                        if flag == 20
                           error('Unexpected problem when initializing data!!')                                 
                            
                        end
                        % map whole session data:
                        if ~exist(fullPath,'file')
                            fid = fopen(fullPath, 'w');
                            fwrite(fid, I, 'single');
                            fclose(fid);clear fid;
                        end
                        memMap.('whole').([Priv.prefixFreqBand ...
                            num2str(m)]).([Priv.prefixSession ...
                            num2str(s)]).(Priv.simM{3}) = ...
                            memmapfile(fullPath,'format',...
                            {'single',[Priv.dataSize(s,1:3) Priv.nrSubjects]...
                            ,'xyzs'},'Writable',true);
                        % map windowed data:
                        %                    if ~exist([fullPath_win])
                        %                      fid = fopen(fullPath_win, 'w');
                        %                      fwrite(fid, Iw, 'single');
                        %                      fclose(fid);clear fid;
                        %                    end
                        %                    memMap.('win').([Priv.prefixFreqBand ...
                        %                    num2str(m)]).([Priv.prefixSession ...
                        %                    num2str(s)]).(Priv.simM{3}) = ...
                        %                    memmapfile(fullPath_win,'format',...
                        %                    {'single',[Priv.dataSize(s,1:3) Priv.nrSubjects]...
                        %                    ,'xyzs'},'Writable',logical(1));
                        flag = 0;
                    catch err
                        if flag == 20
                           error('Unexpected problem when initializing data!!')                                 
                        end
                        disp(err.message)
                        flag = flag + 1;
                    end
                end
            end
        end
        
        
        clear I;
        load([Pub.dataDestination 'memMaps'])
        memMaps = setfield(memMaps,Priv.withinMapName,memMap);
        save([Pub.dataDestination 'memMaps'],'memMaps')
        disp(' ')
        
    case 'syncCurves'
        LL = 12;
        flagS = false;
        disp('Mapping synch. curves')
        for s = 1:Priv.nrSessions
            if Priv.nrTimeIntervals == 0
                memMap.([Priv.prefixSession num2str(s)]) = [];
            else
                I = zeros([Priv.nrTimeIntervals(s),LL+2,...
                    length(Priv.simM),length(Priv.brainAtlases)/2,length(Priv.brainRegions)]);
                %size(I)
                for m = 0:Priv.maxScale + 1
                    disp(['Mapping session ' num2str(s) ', band ' num2str(m)])
                    % map into memory:
                    fullPath = [Priv.resultsDestination Priv.prefixSyncResults Priv.prefixSession ...
                        num2str(s) Priv.prefixFreqBand num2str(m) '.bin'];
                    if exist(fullPath,'file') == 2
                        disp('Data mapped already...');continue
                    end
                    fid = fopen(fullPath, 'w');
                    % write data to a binary file:
                    fwrite(fid, I, 'single');
                    fclose(fid);clear fid;
                    DimsT = size(I);
                    % perform dynamic memory mapping:
                    memMap.([Priv.prefixSession num2str(s)]).([Priv.prefixFreqBand num2str(m)])= ...
                        memmapfile(fullPath,'format',{'single',DimsT(1:end-1) ,'tcsa'},'Writable',true);
                    flagS = true;
                end
            end
        end
        if flagS
            clear I
            load([Pub.dataDestination 'memMaps'])
            memMaps = setfield(memMaps,Priv.synchMapName,memMap);
            save([Pub.dataDestination 'memMaps'],'memMaps')
            disp(' ')
        end
    case 'phaseSyncCurves'
        disp('Mapping phase synch. curves')
        for s = 1:Priv.nrSessions
            I = zeros([Pub.dataSize(s,4),2+2*(Priv.maxScale+2),length(Priv.brainAtlases)/2,length(Priv.brainRegions)]);
            for m = 0:Priv.maxScale + 1
                disp(['Mapping session ' num2str(s) ', band ' num2str(m)])
                % map into memory:
                fullPath = [Priv.phaseDifDestination Priv.prefixPhaseSyncResults Priv.prefixSession ...
                    num2str(s) Priv.prefixFreqBand num2str(m) '.bin'];
                if exist(fullPath,'file') == 2
                    if s == Priv.nrSessions && m == (Priv.maxScale + 1)
                        disp('ISP curve already initialized...');return
                    else
                        continue
                    end
                end
                fid = fopen(fullPath, 'w');
                % write data to a binary file:
                fwrite(fid, I, 'single');
                fclose(fid);clear fid;
                DimsT = size(I);
                % perform dynamic memory mapping:
                memMap.([Priv.prefixSession num2str(s)]).([Priv.prefixFreqBand num2str(m)])= ...
                    memmapfile(fullPath,'format',{'single',DimsT(1:end-1) ,'tca'},'Writable',true);
            end
        end
        clear I
        load([Pub.dataDestination 'memMaps'])
        memMaps = setfield(memMaps,Priv.phaseSynchMapName,memMap);
        save([Pub.dataDestination 'memMaps'],'memMaps')
        disp(' ')
        
        
        
        
        
        
        
    case 'PF'
        nrFreqComps = ((Priv.maxScale+2)^2-(Priv.maxScale+2))/2;
        disp('Mapping sum ZPF results:')
        I = zeros([Priv.dataSize(1,1:3),7],'single');
        for k = 1:length(Priv.simM)
            if (Pub.corOn && strcmp(Priv.simM{k},'cor'))
                for s = 1:Priv.nrSessions
                    %                Iw = single( zeros([Priv.dataSize(1,1:3),8,Priv.nrTimeIntervals(s)]) );
                    for m = 1:nrFreqComps
                        fullPath = [Priv.PFDestination Priv.prefixPF '_' ...
                            Priv.simM{k} '_' Priv.prefixSession num2str(s) '_' ...
                            Priv.transformType Priv.prefixFreqComp num2str(m)];
                        if exist([fullPath '.bin'],'file') == 2
                            if s == Priv.nrSessions && m == nrFreqComps
                                disp('Frequency comparison maps already initialized...');return
                            else
                                continue
                            end
                        end
                        flag = 1;
                        while flag ~= 0
                            try
                                if flag == 20
                                    error('Problems with writing data, memory mapping quitted!')                                    
                                end
                                if ~exist([fullPath '.bin'],'file')
                                    % map whole session data:
                                    fid = fopen([fullPath '.bin'], 'w');
                                    fwrite(fid, I, 'single');
                                    fclose(fid);clear fid;
                                end
                                memMap.('whole').([Priv.prefixSession ...
                                    num2str(s)]).(Priv.simM{k}).([Priv.prefixFreqComp num2str(m)]) = ...
                                    memmapfile([fullPath '.bin'],'format',...
                                    {'single',[Priv.dataSize(s,1:3),7],'xyzc'},'Writable',true);
                                disp(['Session: ' num2str(s) ', freq.band comp: ' num2str(m)])
                                % map windowed data:
                                %                      disp(['Mapping Pearson-Filon matrices ' ...
                                %                      ' session: ' num2str(s) ', freq.band comp: ' num2str(m)])
                                %                      fullPath_win = [fullPath '_win'];
                                %                      if ~exist([fullPath_win '.bin'])
                                %                        fid = fopen([fullPath_win '.bin'], 'w');
                                %                        fwrite(fid, Iw, 'single');
                                %                        fclose(fid);clear fid;
                                %                      end
                                %                      memMap.('win').([Priv.prefixSession ...
                                %                      num2str(s)]).(Priv.simM{k}).([Priv.prefixFreqComp ...
                                %                      num2str(m)]) = ...
                                %                      memmapfile([fullPath_win '.bin'],'format',...
                                %                      {'single',[Priv.dataSize(s,1:3),8]...
                                %                      ,'xyzc'},'Writable',logical(1));
                                flag = 0;
                            catch err
                                if flag == 20
                                    error('Problems with writing data, memory mapping quitted!!')
                                end
                                disp(err.message)
                                flag = flag + 1;
                            end
                        end
                        
                    end
                end
            end
        end
        
        clear I;
        load([Pub.dataDestination 'memMaps'])
        memMaps = setfield(memMaps,Priv.PFMapName,memMap);
        save([Pub.dataDestination 'memMaps'],'memMaps')
        disp(' ')
        
    case 'PFmats'
        nrFreqComps = ((Priv.maxScale+2)^2-(Priv.maxScale+2))/2;
        disp('Mapping sum ZPF matrices:')
        I = zeros([Priv.dataSize(1,1:3),(Priv.nrSubjects^2-Priv.nrSubjects)/2],'single');
        for k = 1:length(Priv.simM)
            if (Pub.corOn && strcmp(Priv.simM{k},'cor'))
                for s = 1:Priv.nrSessions
                    %                Iw = single( zeros([Priv.dataSize(1,1:3),2,Priv.nrTimeIntervals(s)]) );
                    for m = 1:nrFreqComps
                        fullPath = [Priv.PFDestination Priv.prefixPFMat '_' ...
                            Priv.simM{k} '_' Priv.prefixSession num2str(s) '_' ...
                            Priv.transformType Priv.prefixFreqComp num2str(m)];
                        if exist([fullPath '.bin'],'file') == 2
                            if m == nrFreqComps && s == Priv.nrSessions
                                disp('Frequency comparison maps already initialized...');return
                            else
                                continue
                            end
                        end
                        
                        flag = 1;
                        while flag ~= 0
                            try
                                if flag == 20
                                    error('Problems with writing data, memory mapping quitted!')
                                end
                                if ~exist([fullPath '.bin'],'file')
                                    % map whole session data:
                                    fid = fopen([fullPath '.bin'], 'w');
                                    fwrite(fid, I, 'single');
                                    fclose(fid);clear fid;
                                end
                                memMap.('whole').([Priv.prefixSession ...
                                    num2str(s)]).(Priv.simM{k}).([Priv.prefixFreqComp num2str(m)]) = ...
                                    memmapfile([fullPath '.bin'],'format',...
                                    {'single',[Priv.dataSize(s,1:3),(Priv.nrSubjects^2-Priv.nrSubjects)/2]...
                                    ,'xyzc'},'Writable',true);
                                disp(['Session: ' num2str(s) ', freq.band comp: ' num2str(m)])
                                % map windowed data:
                                %                     disp(['Mapping Pearson-Filon matrices ' ...
                                %                     ' session: ' num2str(s) ', freq.band comp: ' num2str(m)])
                                %                     fullPath_win = [fullPath '_win'];
                                %                     if ~exist([fullPath_win '.bin'])
                                %                       fid = fopen([fullPath_win '.bin'], 'w');
                                %                       fwrite(fid, Iw, 'single');
                                %                       fclose(fid);clear fid;
                                %                     end
                                %                      memMap.('win').([Priv.prefixSession ...
                                %                      num2str(s)]).(Priv.simM{k}).([Priv.prefixFreqComp ...
                                %                      num2str(m)]) = ...
                                %                      memmapfile([fullPath_win '.bin'],'format',...
                                %                      {'single',[Priv.dataSize(s,1:3),2]...
                                %                      ,'xyzc'},'Writable',logical(1));
                                flag = 0;
                            catch err
                                if flag == 20
                                    error('Problems with writing data, memory mapping quitted!!')
                                end
                                disp(err.message)
                                flag = flag + 1;
                            end
                        end
                        
                    end
                end
            end
        end
        
        clear I;
        load([Pub.dataDestination 'memMaps'])
        memMaps = setfield(memMaps,Priv.PFmatMapName,memMap);
        save([Pub.dataDestination 'memMaps'],'memMaps')
        disp(' ')
        
    case 'PFSession'
        nrSessComps = ((Priv.nrSessions)^2-(Priv.nrSessions))/2;
        disp('Mapping sum ZPF maps across sessions:')
        I = single( zeros([Priv.dataSize(1,1:3),7]) );
        for k = 1:length(Priv.simM)
            if (Pub.corOn && strcmp(Priv.simM{k},'cor'))
                for s = 0:Pub.nrFreqBands
                    for m = 1:nrSessComps
                        fullPath = [Priv.PFsessionDestination Priv.prefixPF '_' ...
                            Priv.simM{k} '_' Priv.prefixFreqBand num2str(s) '_' ...
                            Priv.transformType Priv.prefixSessComp num2str(m)];
                        if exist([fullPath '.bin']) == 2
                            if s == Pub.nrFreqBands && m == nrSessComps
                                disp('Session comparison maps already initialized...');return
                            else
                                continue
                            end
                        end
                        flag = 1;
                        while flag ~= 0
                            try
                                if flag == 20
                                    error
                                    return
                                end
                                if ~exist([fullPath '.bin'])
                                    % map whole session data:
                                    fid = fopen([fullPath '.bin'], 'w');
                                    fwrite(fid, I, 'single');
                                    fclose(fid);clear fid;
                                end
                                memMap.('whole').([Priv.prefixFreqBand num2str(s)]).(Priv.simM{k}).([...
                                    Priv.prefixSessComp num2str(m)]) = ...
                                    memmapfile([fullPath '.bin'],'format',...
                                    {'single',[Priv.dataSize(1,1:3),7],'xyzc'},'Writable',logical(1));
                                
                                disp(['Freq. band: ' num2str(s) ', session comp: ' num2str(m)])
                                % map windowed data:
                                %                      disp(['Mapping Pearson-Filon matrices ' ...
                                %                      ' session: ' num2str(s) ', freq.band comp: ' num2str(m)])
                                %                      fullPath_win = [fullPath '_win'];
                                %                      if ~exist([fullPath_win '.bin'])
                                %                        fid = fopen([fullPath_win '.bin'], 'w');
                                %                        fwrite(fid, Iw, 'single');
                                %                        fclose(fid);clear fid;
                                %                      end
                                %                      memMap.('win').([Priv.prefixSession ...
                                %                      num2str(s)]).(Priv.simM{k}).([Priv.prefixFreqComp ...
                                %                      num2str(m)]) = ...
                                %                      memmapfile([fullPath_win '.bin'],'format',...
                                %                      {'single',[Priv.dataSize(s,1:3),8]...
                                %                      ,'xyzc'},'Writable',logical(1));
                                flag = 0;
                            catch
                                if flag == 20
                                    error('Problems with writing data, memory mapping quitted!!')
                                    return
                                end
                                disp(lasterr)
                                flag = flag + 1;
                            end
                        end
                        
                    end
                end
            end
        end
        
        clear I;
        load([Pub.dataDestination 'memMaps'])
        memMaps = setfield(memMaps,Priv.PFMapSessionName,memMap);
        save([Pub.dataDestination 'memMaps'],'memMaps')
        disp(' ')
        
    case 'PFmatsSession'                 
        nrSessComps = ((Priv.nrSessions)^2-(Priv.nrSessions))/2;
        disp('Mapping sum ZPF matrices:')
        I = single( zeros([Priv.dataSize(1,1:3),(Priv.nrSubjects^2-Priv.nrSubjects)/2]) );
        for k = 1:length(Priv.simM)
            if (Pub.corOn && strcmp(Priv.simM{k},'cor'))
                %                Iw = single( zeros([Priv.dataSize(1,1:3),2,Priv.nrTimeIntervals(s)]) );
                for s = 0:Pub.nrFreqBands
                    for m = 1:nrSessComps
                        fullPath = [Priv.PFsessionDestination Priv.prefixPFMat '_' ...
                            Priv.simM{k} '_' Priv.prefixFreqBand num2str(s) '_' ...
                            Priv.transformType Priv.prefixSessComp num2str(m)];
                        if exist([fullPath '.bin']) == 2
                            if s == Pub.nrFreqBands && m == nrSessComps
                                disp('Session comparison maps already initialized...');
                                return
                            else                          
                                continue
                            end
                        end
                                                                                                
                        flag = 1;
                        while flag ~= 0
                            try
                                if flag == 20
                                    error
                                    return
                                end
                                if ~exist([fullPath '.bin'])
                                    % map whole session data:
                                    fid = fopen([fullPath '.bin'], 'w');
                                    fwrite(fid, I, 'single');
                                    fclose(fid);clear fid;
                                end
                                memMap.('whole').([Priv.prefixFreqBand ...
                                    num2str(s)]).(Priv.simM{k}).([Priv.prefixSessComp num2str(m)]) = ...
                                    memmapfile([fullPath '.bin'],'format',...
                                    {'single',[Priv.dataSize(1,1:3),(Priv.nrSubjects^2-Priv.nrSubjects)/2]...
                                    ,'xyzc'},'Writable',logical(1));
                                disp(['Freq. band: ' num2str(s) ', session comp: ' num2str(m)])
                                
                                
                                % map windowed data:
                                %                     disp(['Mapping Pearson-Filon matrices ' ...
                                %                     ' session: ' num2str(s) ', freq.band comp: ' num2str(m)])
                                %                     fullPath_win = [fullPath '_win'];
                                %                     if ~exist([fullPath_win '.bin'])
                                %                       fid = fopen([fullPath_win '.bin'], 'w');
                                %                       fwrite(fid, Iw, 'single');
                                %                       fclose(fid);clear fid;
                                %                     end
                                %                      memMap.('win').([Priv.prefixSession ...
                                %                      num2str(s)]).(Priv.simM{k}).([Priv.prefixFreqComp ...
                                %                      num2str(m)]) = ...
                                %                      memmapfile([fullPath_win '.bin'],'format',...
                                %                      {'single',[Priv.dataSize(s,1:3),2]...
                                %                      ,'xyzc'},'Writable',logical(1));
                                flag = 0;
                            catch
                                if flag == 20
                                    error('Problems with writing data, memory mapping quitted!!')
                                    return
                                end
                                disp(lasterr)
                                flag = flag + 1;
                            end
                        end
                        
                    end
                end
            end
        end
        
        clear I;
        load([Pub.dataDestination 'memMaps'])
        memMaps = setfield(memMaps,Priv.PFmatMapSessionName,memMap);
        save([Pub.dataDestination 'memMaps'],'memMaps')
        disp(' ')
        
        
        
        
        
    case 'cormats'
        disp('Mapping ISC matrices:')
        I = zeros([Priv.dataSize(1,1:3),(Priv.nrSubjects^2-Priv.nrSubjects)/2],'single');
        for k = 1:length(Priv.simM)
            if (Pub.corOn && strcmp(Priv.simM{k},'cor'))
                for s = 1:Priv.nrSessions
                    for m = 0:Priv.maxScale + 1
                        
                        disp(['Mapping stats ' ...
                            ' session ' num2str(s) ', band ' num2str(m)])
                        
                        fullPath = [Priv.statsDestination ...
                            Priv.prefixCorMat '_' ...
                            Priv.simM{k} '_' Priv.prefixFreqBand ...
                            num2str(m) '_' Priv.prefixSession ...
                            num2str(s) '_' Priv.transformType];
                        if exist([fullPath '.bin'],'file') == 2
                            disp('Data already exist, quit memory mapping...');return
                        end
                        flag = 1;
                        while flag ~= 0
                            try
                                if flag == 20
                                    error('Problems with writing data, memory mapping quitted!')
                                end
                                % map whole session data:
                                if ~exist([fullPath '.bin'],'file')
                                    fid = fopen([fullPath '.bin'], 'w');
                                    fwrite(fid, I, 'single');
                                    fclose(fid);clear fid;
                                end
                                memMap.('whole').([Priv.prefixFreqBand ...
                                    num2str(m)]).([Priv.prefixSession ...
                                    num2str(s)]).(Priv.simM{k}) = ...
                                    memmapfile([fullPath '.bin'],'format',...
                                    {'single',[Priv.dataSize(s,1:3),...
                                    (Priv.nrSubjects^2-Priv.nrSubjects)/2]...
                                    ,'xyzc'},'Writable',true);
                                % map windowed data:
                                for tv = 1:Priv.nrTimeIntervals(s)
                                    fullPath_win = [fullPath '_win' num2str(tv)];
                                    if exist([fullPath_win '.bin'],'file') == 2
                                        disp('Data already exist, quit memory mapping...');return
                                    end
                                    if ~exist([fullPath_win '.bin'],'file')
                                        fid = fopen([fullPath_win '.bin'], 'w');
                                        fwrite(fid, I, 'single');
                                        fclose(fid);clear fid;
                                    end
                                    memMap.('win').([Priv.prefixFreqBand ...
                                        num2str(m)]).([Priv.prefixSession ...
                                        num2str(s)]).(Priv.simM{k}).([Priv.prefixTimeVal num2str(tv)]) = ...
                                        memmapfile([fullPath_win '.bin'],'format',...
                                        {'single',[Priv.dataSize(s,1:3),...
                                        (Priv.nrSubjects^2-Priv.nrSubjects)/2]...
                                        ,'xyzc'},'Writable',true);
                                    flag = 0;
                                end
                                flag = 0;
                            catch err
                                if flag == 20
                                    error('Problems with writing data, memory mapping quitted!!')
                                end
                                disp(err.message)
                                flag = flag + 1;
                            end
                            
                        end
                    end
                end
            end
        end
        clear I;
        load([Pub.dataDestination 'memMaps'])
        memMaps = setfield(memMaps,Priv.cormatMapName,memMap);
        save([Pub.dataDestination 'memMaps'],'memMaps')
        disp(' ')
        
    case 'phase'
        disp('Mapping phase synchronization data:')
        for s = 1:Priv.nrSessions
            
            for m = 0:Priv.maxScale + 1
                fullPath = [Priv.phaseDifDestination ...
                    Priv.prefixPhaseDif '_' Priv.prefixFreqBand ...
                    num2str(m) '_' Priv.prefixSession ...
                    num2str(s) '_' Priv.transformType '.bin'];
                if exist(fullPath,'file') == 2
                    disp('Data already exist, quit memory mapping...');return
                end
                I = zeros(Priv.dataSize(s,:),'single');
                I = permute(I,[4 2 3 1]);
                flag = 1;
                while flag ~= 0
                    try
                        if flag == 20
                            error('Problems with writing data, memory mapping quitted!')
                            
                        end
                        if ~exist([fullPath '.bin'],'file')
                            fid = fopen(fullPath, 'w');
                            % write data to a binary file:
                            fwrite(fid, I, 'single');
                            fclose(fid);clear fid;
                        end
                        % perform dynamic memory mapping:
                        memMap.([Priv.prefixSession ...
                            num2str(s)]).([Priv.prefixFreqBand ...
                            num2str(m)]) = memmapfile(fullPath, ...
                            'format',{'single',[Priv.dataSize(...
                            s,4) Priv.dataSize(s,2:3)]...
                            ,'tyz'},'Writable',true);
                        disp(['Session ' num2str(s)  ', Frequency band ' num2str(m)])
                        flag = 0;
                    catch err
                        if flag == 20
                            error('Problems with writing data, memory mapping quitted!!')
                            
                        end
                        disp(err.message)
                        flag = flag + 1;
                    end
                end
            end
        end
        clear I;
        load([Pub.dataDestination 'memMaps'])
        memMaps = setfield(memMaps,Priv.phaseMapName,memMap);
        save([Pub.dataDestination 'memMaps'],'memMaps')
        disp(' ')
        
    case 'stats'
        flagStats = [false];
        disp('Mapping other ISC based maps:')
        I = zeros([Priv.dataSize(1,1:3),5],'single' );
        for k = 1:length(Priv.simM)
            if (Pub.corOn && strcmp(Priv.simM{k},'cor'))
                for s = 1:Priv.nrSessions
                    Iw = zeros([Priv.dataSize(1,1:3),5,Priv.nrTimeIntervals(s)],'single');
                    for m = 0:Priv.maxScale + 1
                        disp(['Mapping stats ' ...
                            ' session ' num2str(s) ', band ' num2str(m)])
                            fullPath = [Priv.statsDestination ...
                            Priv.prefixTMap '_' ...
                            Priv.simM{k} '_' Priv.prefixFreqBand ...
                            num2str(m) '_' Priv.prefixSession ...
                            num2str(s) '_' Priv.transformType];
                            if Pub.winOn
                                fullPath_win = [fullPath '_win'];
                            end
                            if Pub.winOn
                                if exist([fullPath '.bin'],'file') == 2 && exist([fullPath_win '.bin'],'file') == 2                                    
                                    if s == Priv.nrSessions && (m == Priv.maxScale+1)
                                        flagStats(1) = true;
                                    else
                                        continue
                                    end
                                end
                            else
                                if exist([fullPath '.bin'],'file') == 2 
                                    if s == Priv.nrSessions && (m == Priv.maxScale+1)
                                        flagStats(1) = true;
                                    else
                                        continue
                                    end
                                end
                            end                      
                        if isequal(flagStats,true)
                           disp('Other ISC maps already intialized...')
                           return
                        end
                        
                        flag = 1;
                        while flag ~= 0
                            try
                                if flag == 20
                                    error('Unexpected problem when initializing data!!')                                   
                                end
                                % map whole session data:
                                if ~exist([fullPath '.bin'],'file')
                                    fid = fopen([fullPath '.bin'], 'w');
                                    fwrite(fid, I, 'single');
                                    fclose(fid);clear fid;
                                end
                                memMap.('whole').([Priv.prefixFreqBand ...
                                    num2str(m)]).([Priv.prefixSession ...
                                    num2str(s)]).(Priv.simM{k}) = ...
                                    memmapfile([fullPath '.bin'],'format',...
                                    {'single',[Priv.dataSize(s,1:3),5]...
                                    ,'xyz'},'Writable',true);
                                % map windowed data:
                                if Pub.winOn
                                    if ~exist([fullPath_win '.bin'],'file')
                                        fid = fopen([fullPath_win '.bin'], 'w');
                                        fwrite(fid, Iw, 'single');
                                        fclose(fid);clear fid;
                                    end
                                    memMap.('win').([Priv.prefixFreqBand ...
                                        num2str(m)]).([Priv.prefixSession ...
                                        num2str(s)]).(Priv.simM{k}) = ...
                                        memmapfile([fullPath_win '.bin'],'format',...
                                        {'single',[Priv.dataSize(s,1:3),5]...
                                        ,'xyz'},'Writable',true);
                                end
                                flag = 0;
                                
                            catch err
                                if flag == 20
                                   error('Unexpected problem when initializing data!!') 
                                else
                                    disp(err.message)
                                    flag = flag + 1;
                                end
                            end
                            
                        end
                    end
                end
            end
        end
        clear I;
        load([Pub.dataDestination 'memMaps'])
        memMaps = setfield(memMaps,Priv.statMapName,memMap);
        save([Pub.dataDestination 'memMaps'],'memMaps')
        disp(' ')
        
    otherwise
        error('Unknown data specification!!')
        
end

function ParamsNew = checkUpdates(Params,ParamsNew)

PubNew = ParamsNew.PublicParams;
PrivNew = ParamsNew.PrivateParams;
Pub = Params.PublicParams;
Priv = Params.PrivateParams;


if ( ( PubNew.windowSize ~= Pub.windowSize ) || ( PubNew.windowStep ~= Pub.windowStep) )
    PrivNew.filterUpdate = false;
    disp('Changed time-window parameters. New results will overwrite the existing ones in the directory.')
end

if ( ( PubNew.ssiOn ~= Pub.ssiOn ) || ( PubNew.nmiOn ~= Pub.nmiOn ) || ...
        ( PubNew.corOn ~= Pub.corOn ) || ( PubNew.kenOn ~= Pub.kenOn ) )
    PrivNew.filterUpdate = false;
    disp('Changed parameters found for similarity maps, specified maps will overwrite the existing ones in the directory.')
end

ParamsNew.PrivateParams = PrivNew;
ParamsNew.PublicParams = PubNew;
