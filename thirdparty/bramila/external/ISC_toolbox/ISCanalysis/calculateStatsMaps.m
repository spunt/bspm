function calculateStatsMaps(Params,nrBand,nrSession)

% This function calculates time-varying and across-session ISC maps
% maps based on the following statistic:
% 1. The mean correlation coefficient after applying Fisher's
% z-tranformation to subject pairwise correlation coefficients.
% 2. The corresponding standard deviation.
% 3. The corresponding lower-quartile.
% 4. The corresponding median.
% 5. The corresponding upper-quartile.
%
% Results can be accessed through memMaps.statMap -field.
% Inputs:
% Params - struct containing all necessary parameters
% nrBand - frequency subband index (note: 0 refers to full band)
% nrSession - session index
%
% Example:
% load analysisParameters
% calculateStatsMaps(analysisParameters,4,1);
% Then access maps of 10th time-interval:
% load memMaps
% statMaps = memMaps.statMap.win.band4.Session1.cor.Data(10).xyz;
% meanFisMap = statMaps(:,:,1);
% stdFisMap = statMaps(:,:,2);
% lowQuartMap = statMaps(:,:,3);
% medianMap = statMaps(:,:,4);
% upQuartMap = statMaps(:,:,5);
%
% See also:
% ISCANALYSIS
% RUNANALYSIS
% INITPARAMS
% CALCULATESIMILARITYMAPS
% CALCULATECORMATS

% Last modified 5.8.2013 by Juha Pajula
% Tampere University of Technology
% Department of Signal Processing
% e-mail: juha.pajula@tut.fi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

showTime(1);

Pub = Params.PublicParams;
Priv = Params.PrivateParams;

if ~Pub.calcStats
    disp('Additional ISC maps not selected...')
    return
end
if Pub.corOn == 0
    disp('Additional ISC maps can be computed only when basic ISC analysis is selected...')
    return
end

nrSubjectPairs = Priv.nrSubjects*(Priv.nrSubjects-1)/2;
if nrSubjectPairs < 4
    warning('Less than four subject pairs found -> calculation of upper and lower quartile skipped...')
end


% load required memory map pointers:
load([Pub.dataDestination 'memMaps'])

if Pub.nrFreqBands > 0
    mMapFilt = memMaps.(Priv.filtMapName).([Priv.prefixSession num2str(nrSession)]);
end

mMapOrig = memMaps.(Priv.origMapName).([Priv.prefixSession num2str(nrSession)]);
mMapStat = memMaps.(Priv.statMapName);
clear memMaps
try
    if ~isempty(mMapStat)
        if isfield(mMapStat,'whole')
            if isempty(mMapStat.whole) == false
                if mMapStat.whole.([Priv.prefixFreqBand num2str(nrBand)]).(...
                        [Priv.prefixSession num2str(nrSession)]).cor.Writable == false
                    disp('Additional ISC maps computed already...')
                    return
                end
            end
        end
        if isfield(mMapStat,'win')
            if isempty(mMapStat.win) == false
                if mMapStat.win.([Priv.prefixFreqBand num2str(nrBand)]).(...
                        [Priv.prefixSession num2str(nrSession)]).cor.Writable == false
                    disp('Additional ISC maps computed already...')
                    return
                end
            end
        end
    else
        disp('Memory map pointers not available, calculation canceled...')
        return
    end
catch err
    error('Cannot update stat maps possibly due to changes in analysis parameters. To update the data, delete bin-files in stats-folder and try analysis again!')
end


[~,~,en] = computer;

cDat = zeros([Priv.dataSize(nrSession,[4 2 3]), Priv.nrSubjects]);
% load brain mask:
if strcmp(Pub.fileFormatSubj,'nii')
    bmask = load_nii(Priv.brainMask);
    bmask = single(bmask.img);
elseif strcmp(Pub.fileFormatSubj,'mat')
    bmask = load(Priv.brainMask);
    fiel = fields(bmask);
    bmask = bmask.(fiel{1});
    bmask = single(bmask);
else
    error('Mask must be mat- or nii-file!')
end
bmask = logical(bmask);

INDS = find(triu(ones(Priv.nrSubjects,Priv.nrSubjects),1));
iter = 0;

if Priv.nrTimeIntervals(nrSession) == 0
    nrFrames = 1;
else
    nrFrames = Priv.nrTimeIntervals(nrSession);
end

tic
for xx = 1:Priv.dataSize(nrSession,1)
    disp(['x:' num2str(xx) '/' num2str(Priv.dataSize(nrSession,1))])
    % process only non-zero slices:
    if sum(sum(squeeze(bmask(xx,:,:)))) > 0
        cDat = zeros(Priv.dataSize(nrSession,4),Priv.dataSize(nrSession,2),...
            Priv.dataSize(nrSession,3),Priv.nrSubjects);
        
        % get mapped source data of the subjects:
        for k = 1:Priv.nrSubjects
            if nrBand == 0 % load full-band data
                if ~strcmp(Priv.computerInfo.endian,en)
                    cDat(:,:,:,k) = swapbytes(mMapOrig.([Priv.prefixSubject num2str(k)]).Data(xx).tyz);
                else
                    cDat(:,:,:,k) = mMapOrig.([Priv.prefixSubject num2str(k)]).Data(xx).tyz;
                end
            else % load sub-band data
                if ~strcmp(Priv.computerInfo.endian,en)
                    cDat(:,:,:,k) = swapbytes(mMapFilt.([Priv.prefixSubjectFilt ...
                        num2str(k)]).([Priv.prefixFreqBand num2str(nrBand)]).Data(xx).tyz);
                else
                    cDat(:,:,:,k) = mMapFilt.([Priv.prefixSubjectFilt ...
                        num2str(k)]).([Priv.prefixFreqBand num2str(nrBand)]).Data(xx).tyz(:,:,:);
                end
            end
        end
        
        statData = zeros([Priv.dataSize(nrSession,2:3),5]);
        statDatawin = zeros([Priv.dataSize(nrSession,2:3),5,Priv.nrTimeIntervals(nrSession)]);
        
        for yy = 1:Priv.dataSize(nrSession,2)
            for zz = 1:Priv.dataSize(nrSession,3)
                if bmask(xx,yy,zz)
                    iter = iter + 1;
                    if mod(iter,10000) == 0
                        disp(['iter: ' num2str(iter) '/' num2str(length(find(bmask)))])
                        toc
                        tic
                    end
                    % obtain each subject's time series:
                    ts = squeeze(cDat(:,yy,zz,:));
                    % ikkunointi:
                    for wfr = 1:nrFrames
                        % calculate across-session statistical maps:
                        if wfr == 1
                            %N = Priv.dataSize(nrSession,4);
                            if Pub.corOn
                                % correlation coefficient calculation (ref: Matlab corrcoef.m)
                                [n1,m1] = size(ts);
                                xc = ts - repmat(sum(ts)/n1,n1,1);  % Remove mean
                                c1 = (xc' * xc) / (n1-1); % calculate inner products
                                d1 = sqrt(diag(c1)); % sqrt first to avoid under/overflow
                                dd = d1*d1';
                                dd(1:m1+1:end) = diag(c1); % remove roundoff on diag
                                r1 = c1 ./ dd;
                                r1 = r1(INDS);
                                if ~isnan(sum(r1))
                                    % z-transform pairwise correlation values:
                                    fishZ = 0.5*(log((1+r1)./(1-r1)));
                                    % calculate t-map:
                                    sDev = std(fishZ);
                                    if nrSubjectPairs == 1
                                        statData(yy,zz,1) = fishZ;
                                    else
                                        statData(yy,zz,1) = sqrt(length(r1))*mean(fishZ)/sDev;
                                    end
                                    statData(yy,zz,2) = sDev;
                                    % calculate 25%, 50% and 75% percent quartiles:
                                    
                                    % median:
                                    r1S = sort(r1);
                                    nCompare = numel(r1S);
                                    if nCompare <= 1
                                        half = 1;
                                    else
                                        half = floor(nCompare/2);
                                    end
                                    if nCompare > 1
                                        V1 = r1S(half);
                                        V2 = r1S(half+1);
                                        Q2 = V2;
                                        % Average if even number of elements avoiding overflows:
                                        if 2*half == nCompare
                                            Q2 = V1 + (V2-V1)/2;
                                            Id = (sign(V1) ~= sign(V2)) | isinf(V1) | isinf(V2);
                                            Q2(Id) = (V1(Id)+V2(Id))/2;
                                        end
                                    else
                                        Q2 = fishZ;
                                    end
                                    statData(yy,zz,4) = single(Q2);
                                    
                                    % 25% quartile:
                                    if size(r1,1) >= 4
                                        r1S2 = r1S(r1S < Q2);
                                        nCompare = numel(r1S2);
                                        if nCompare <= 1
                                            half = 1;
                                            Q1 = Q2;
                                        else
                                            half = floor(nCompare/2);
                                            V1 = r1S2(half);
                                            V2 = r1S2(half+1);
                                            Q1 = V2;
                                        end
                                        % Average if even number of elements avoiding overflows:
                                        if 2*half == nCompare
                                            Q1 = V1 + (V2-V1)/2;
                                            Id = (sign(V1) ~= sign(V2)) | isinf(V1) | isinf(V2);
                                            Q1(Id) = (V1(Id)+V2(Id))/2;
                                        end
                                        statData(yy,zz,3) = single(Q1);
                                        
                                        % 75% quartile:
                                        r1S3 = r1S(r1S > Q2);
                                        nCompare = numel(r1S3);
                                        if nCompare <= 1
                                            half = 1;
                                            Q3 = Q2;
                                        else
                                            half = floor(nCompare/2);
                                            V1 = r1S3(half);
                                            V2 = r1S3(half+1);
                                            Q3 = V2;
                                        end
                                        % Average if even number of elements avoiding overflows:
                                        if 2*half == nCompare
                                            Q3 = V1 + (V2-V1)/2;
                                            Id = (sign(V1) ~= sign(V2)) | isinf(V1) | isinf(V2);
                                            Q3(Id) = (V1(Id)+V2(Id))/2;
                                        end
                                        statData(yy,zz,5) = single(Q3);
                                    end
                                    
                                    
                                    if mod(iter,10000) == 0
                                        disp(['Intersubject correlation value ('...
                                            num2str(xx) ',' num2str(yy) ',' ...
                                            num2str(zz) '): '])
                                        disp(['t-stat, std, Q25, Q50, Q75: ' ...
                                            num2str(statData(yy,zz,:))])
                                    end
                                else
                                    statData(yy,zz,1) = NaN;
                                    statData(yy,zz,2) = NaN;
                                    statData(yy,zz,3) = NaN;
                                    statData(yy,zz,4) = NaN;
                                    statData(yy,zz,5) = NaN;
                                end
                                
                            end
                            
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % calculate statistical maps for time-frames:
                        if Priv.nrTimeIntervals(nrSession) > 0
                            %N = Pub.windowSize;
                            ts_win = ts(Priv.startInds{nrSession}...
                                (wfr):Priv.startInds{nrSession}...
                                (wfr)+Pub.windowSize-1,:);
                            if Pub.corOn
                                % Mean of pairwise correlation (ref: Matlab corrcoef.m):
                                [n1,m1] = size(ts_win);
                                xc = ts_win - repmat(sum(ts_win)/n1,n1,1);  % Remove mean
                                c1 = (xc' * xc) / (n1-1); % calculate inner products
                                d1 = sqrt(diag(c1)); % sqrt first to avoid under/overflow
                                dd = d1*d1';
                                dd(1:m1+1:end) = diag(c1); % remove roundoff on diag
                                r1 = c1 ./ dd;
                                r1 = r1(INDS);
                                r1 = min(r1,1);
                                r1 = max(r1,-1);
                                if ~isnan(sum(r1))
                                    % z-transform pairwise correlation values:
                                    fishZ = 0.5*(log((1+r1)./(1-r1)));
                                    % calculate t-map:
                                    sDev = std(fishZ);
                                    if nrSubjectPairs == 1
                                        statDatawin(yy,zz,1,wfr) = fishZ;
                                    else
                                        statDatawin(yy,zz,1,wfr) = sqrt(length(r1))*mean(fishZ)/sDev;
                                    end
                                    statDatawin(yy,zz,2,wfr) = sDev;
                                    % calculate 25%, 50% and 75% percent quartiles:
                                    % median:
                                    
                                    r1S = sort(r1);
                                    nCompare = numel(r1S);
                                    if nCompare > 1
                                        half = floor(nCompare/2);
                                        V1 = r1S(half);
                                        V2 = r1S(half+1);
                                        Q2 = V2;
                                        % Average if even number of elements avoiding overflows:
                                        if 2*half == nCompare
                                            Q2 = V1 + (V2-V1)/2;
                                            Id = (sign(V1) ~= sign(V2)) | isinf(V1) | isinf(V2);
                                            Q2(Id) = (V1(Id)+V2(Id))/2;
                                        end
                                    else
                                        Q2 = fishZ;
                                    end
                                    
                                    statDatawin(yy,zz,4,wfr) = single(Q2);
                                    % 25% quartile:
                                    if size(r1,1) >= 4
                                        r1S2 = r1S(r1S < Q2);
                                        nCompare = numel(r1S2);
                                        if nCompare <= 1
                                            half = 1;
                                            Q1 = Q2;
                                        else
                                            half = floor(nCompare/2);
                                            V1 = r1S2(half);
                                            V2 = r1S2(half+1);
                                            Q1 = V2;
                                        end
                                        % Average if even number of elements avoiding overflows:
                                        if 2*half == nCompare
                                            Q1 = V1 + (V2-V1)/2;
                                            Id = (sign(V1) ~= sign(V2)) | isinf(V1) | isinf(V2);
                                            Q1(Id) = (V1(Id)+V2(Id))/2;
                                        end
                                        % 75% quartile:
                                        r1S3 = r1S(r1S > Q2);
                                        nCompare = numel(r1S3);
                                        if nCompare <= 1
                                            half = 1;
                                            Q3 = Q2;
                                        else
                                            half = floor(nCompare/2);
                                            V1 = r1S3(half);
                                            V2 = r1S3(half+1);
                                            Q3 = V2;
                                        end
                                        % Average if even number of elements avoiding overflows:
                                        if 2*half == nCompare
                                            Q3 = V1 + (V2-V1)/2;
                                            Id = (sign(V1) ~= sign(V2)) | isinf(V1) | isinf(V2);
                                            Q3(Id) = (V1(Id)+V2(Id))/2;
                                        end
                                        statDatawin(yy,zz,3,wfr) = single(Q1);
                                        statDatawin(yy,zz,5,wfr) = single(Q3);
                                    end
                                else
                                    statDatawin(yy,zz,1,wfr) = NaN;
                                    statDatawin(yy,zz,2,wfr) = NaN;
                                    statDatawin(yy,zz,3,wfr) = NaN;
                                    statDatawin(yy,zz,4,wfr) = NaN;
                                    statDatawin(yy,zz,5,wfr) = NaN;
                                end
                                
                            end
                        end
                    end
                end
            end
        end
        
        % save results:
        %      tic
        if Pub.corOn
            mMapStat.whole.([Priv.prefixFreqBand...
                num2str(nrBand)]).([Priv.prefixSession...
                num2str(nrSession)]).cor.Data.xyz(xx,:,:,:) = statData;
        end
        if Priv.nrTimeIntervals(nrSession) > 0
            for wfr = 1:size(statDatawin,4)
                if Pub.corOn
                    mMapStat.win.([Priv.prefixFreqBand...
                        num2str(nrBand)]).([Priv.prefixSession...
                        num2str(nrSession)]).cor.Data(wfr).xyz(xx,:,:,:) = statDatawin(:,:,:,wfr);
                end
            end
        end
        %      toc
    end
end

if Pub.corOn
    %    lock_name=['calcStatsMaps_',num2str(subjecNr),'_',num2str(sessionNr)];
    lock_name=['calcStatsMaps_',num2str(nrBand),'_',num2str(nrSession)];
    
    if(freeToWrite('check',Pub.dataDestination,lock_name))
        
        load([Pub.dataDestination 'memMaps'])
        memMaps.(Priv.statMapName).whole.([Priv.prefixFreqBand...
            num2str(nrBand)]).([Priv.prefixSession...
            num2str(nrSession)]).cor.Writable = false;
        save([Pub.dataDestination 'memMaps.mat'],'memMaps')
        clear memMaps
        if Priv.nrTimeIntervals(nrSession) > 0
            load([Pub.dataDestination 'memMaps'])
            memMaps.(Priv.statMapName).win.([Priv.prefixFreqBand...
                num2str(nrBand)]).([Priv.prefixSession...
                num2str(nrSession)]).cor.Writable = false;
            save([Pub.dataDestination 'memMaps.mat'],'memMaps')
            clear memMaps
        end
        [~]=freeToWrite('release',Pub.dataDestination,lock_name);
    end
end


showTime(0);
