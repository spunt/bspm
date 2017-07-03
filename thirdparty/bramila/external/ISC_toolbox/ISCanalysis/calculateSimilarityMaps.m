function calculateSimilarityMaps(Params,nrBand,nrSession)

% This function calculates time-varying and across-session
% intersubject similarity maps. Results can be accessed
% through memMaps.resultMap -field.
%
% Inputs:
% Params - struct containing all necessary parameters
% nrBand - frequency subband index (note: 0 refers to full band)
% nrSession - session index
%
% Example:
% load analysisParameters
% calculateSimilarityMaps(analysisParameters,4,1);
% Then access ISC-map of 10th time-interval:
% load memMaps
% ISCmap10 = memMaps.resultMap.win.band4.Session1.cor.Data(10).xyz;
%
% See also:
% INITPARAMS
% RUNANALYSIS


% Last modified 1.11.2013 by Jukka-Pekka Kauppi
% University of Helsinki
% Department of Computer Science
% e-mail: jukka-pekka.kauppi@helsinki.fi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


showTime(1);

Pub = Params.PublicParams;
Priv = Params.PrivateParams;


% check parameters:

if ~Pub.calcStandard
    disp('Basic ISC analysis not selected...')
    return
end

load([Pub.dataDestination 'memMaps'])

if Pub.nrFreqBands > 0
    mMapFilt = memMaps.(Priv.filtMapName).([Priv.prefixSession num2str(nrSession)]);
end

mMapOrig = memMaps.(Priv.origMapName).([Priv.prefixSession num2str(nrSession)]);
mMapResult = memMaps.(Priv.resultMapName);
flag = false*ones(1,8);
%try
if ~isempty(mMapResult)
    if isfield(mMapResult,'whole')
        if isempty(mMapResult.whole) == false
            if isfield(mMapResult.whole,[Priv.prefixFreqBand num2str(nrBand)])
                if isfield(mMapResult.whole.([Priv.prefixFreqBand num2str(nrBand)]),...
                        ([Priv.prefixSession num2str(nrSession)]))
                    if isfield(mMapResult.whole.([Priv.prefixFreqBand num2str(nrBand)])...
                            .([Priv.prefixSession num2str(nrSession)]),'cor')
                        if Pub.corOn
                            if mMapResult.whole.([Priv.prefixFreqBand num2str(nrBand)]).(...
                                    [Priv.prefixSession num2str(nrSession)]).cor.Writable == false
                                disp('ISC maps already computed...')
                                flag(1) = true;
                            end
                        end
                    end
                    if isfield(mMapResult.whole.([Priv.prefixFreqBand num2str(nrBand)])...
                            .([Priv.prefixSession num2str(nrSession)]),'ssi')
                        if Pub.ssiOn
                            if mMapResult.whole.([Priv.prefixFreqBand num2str(nrBand)]).(...
                                    [Priv.prefixSession num2str(nrSession)]).ssi.Writable == false
                                disp('SSI maps already computed...')
                                flag(2) = true;
                            end
                        end
                    end
                    if isfield(mMapResult.whole.([Priv.prefixFreqBand num2str(nrBand)])...
                            .([Priv.prefixSession num2str(nrSession)]),'nmi')
                        if Pub.nmiOn
                            if mMapResult.whole.([Priv.prefixFreqBand num2str(nrBand)]).(...
                                    [Priv.prefixSession num2str(nrSession)]).nmi.Writable == false
                                disp('NMI maps already computed...')
                                flag(3) = true;
                            end
                        end
                    end
                    if isfield(mMapResult.whole.([Priv.prefixFreqBand num2str(nrBand)])...
                            .([Priv.prefixSession num2str(nrSession)]),'ken')
                        if Pub.kenOn
                            if mMapResult.whole.([Priv.prefixFreqBand num2str(nrBand)]).(...
                                    [Priv.prefixSession num2str(nrSession)]).ken.Writable == false
                                disp('Kendall maps already computed...')
                                flag(4) = true;
                            end
                        end
                    end
                end
            end
        end
    end
    if isfield(mMapResult,'win')
        if isempty(mMapResult.win) == false
            if isfield(mMapResult.win,[Priv.prefixFreqBand num2str(nrBand)])
                if isfield(mMapResult.win.([Priv.prefixFreqBand num2str(nrBand)]),...
                        ([Priv.prefixSession num2str(nrSession)]))
                    if isfield(mMapResult.win.([Priv.prefixFreqBand num2str(nrBand)])...
                            .([Priv.prefixSession num2str(nrSession)]),'cor')
                        if Pub.corOn
                            if mMapResult.win.([Priv.prefixFreqBand num2str(nrBand)]).(...
                                    [Priv.prefixSession num2str(nrSession)]).cor.Writable == false
                                disp('Windowed ISC maps already computed...')
                                flag(5) = true;
                            end
                        end
                    end
                    if isfield(mMapResult.win.([Priv.prefixFreqBand num2str(nrBand)])...
                            .([Priv.prefixSession num2str(nrSession)]),'ssi')
                        if Pub.ssiOn
                            if mMapResult.win.([Priv.prefixFreqBand num2str(nrBand)]).(...
                                    [Priv.prefixSession num2str(nrSession)]).ssi.Writable == false
                                disp('Windowed SSI maps already computed...')
                                flag(6) = true;
                            end
                        end
                    end
                    if isfield(mMapResult.win.([Priv.prefixFreqBand num2str(nrBand)])...
                            .([Priv.prefixSession num2str(nrSession)]),'nmi')
                        if Pub.nmiOn
                            if mMapResult.win.([Priv.prefixFreqBand num2str(nrBand)]).(...
                                    [Priv.prefixSession num2str(nrSession)]).nmi.Writable == false
                                disp('Windowed NMI maps already computed...')
                                flag(7) = true;
                            end
                        end
                    end
                    if isfield(mMapResult.win.([Priv.prefixFreqBand num2str(nrBand)])...
                            .([Priv.prefixSession num2str(nrSession)]),'ken')
                        
                        if Pub.kenOn
                            if mMapResult.win.([Priv.prefixFreqBand num2str(nrBand)]).(...
                                    [Priv.prefixSession num2str(nrSession)]).ken.Writable == false
                                disp('Windowed Kendall maps already computed...')
                                flag(8) = true;
                            end
                        end
                    end
                end
            end
        end
    end
else
    disp('Memory pointers not available, calculation canceled...')
    return
end
%catch err
%   error('Cannot update similarity maps possibly due to changes in analysis parameters. To update the data, delete bin-files in results-folder and try analysis again!')
%end

clear memMaps

selectedAnalysis = false*ones(1,8);
if Pub.corOn
    selectedAnalysis(1) = true;
end
if Pub.ssiOn
    selectedAnalysis(2) = true;
end
if Pub.nmiOn
    selectedAnalysis(3) = true;
end
if Pub.kenOn
    selectedAnalysis(4) = true;
end
if Pub.winOn
    selectedAnalysis(5:8) = selectedAnalysis(1:4);
end
if isequal(selectedAnalysis,flag)
    disp('All selected similarity maps computed...')
    return
end

[~,~,en] = computer;

cDat = zeros([Priv.dataSize(nrSession,[4 2 3]), Priv.nrSubjects]);

% load mask:
if strcmp(Pub.fileFormat,'nii')
    bmask = load_nii(Priv.brainMask);
    bmask = single(bmask.img);
elseif strcmp(Pub.fileFormat,'mat')
    bmask = load(Priv.brainMask);
    fiel = fields(bmask);
    bmask = bmask.(fiel{1});
    bmask = single(bmask);
else
    error('Mask must be mat- or nii-file!')
end
bmask = logical(bmask);

INDS = find(triu(ones(Priv.nrSubjects,Priv.nrSubjects),1));
% INDS = [];
% for hh = 1:Priv.nrSubjects-1
%     INDS = [INDS (1+hh*Priv.nrSubjects):(...
%         1+hh*Priv.nrSubjects+hh-1)];
% end
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
        
        % init temp data matrices:
        kenDatawin = zeros([Priv.dataSize(nrSession,2:3),Priv.nrTimeIntervals(nrSession)]);
        corDatawin = zeros([Priv.dataSize(nrSession,2:3),Priv.nrTimeIntervals(nrSession)]);
        ssiDatawin = zeros([Priv.dataSize(nrSession,2:3),Priv.nrTimeIntervals(nrSession)]);
        nmiDatawin = zeros([Priv.dataSize(nrSession,2:3),Priv.nrTimeIntervals(nrSession)]);
        kenData = zeros(Priv.dataSize(nrSession,2:3));
        corData = zeros(Priv.dataSize(nrSession,2:3));
        ssiData = zeros(Priv.dataSize(nrSession,2:3));
        nmiData = zeros(Priv.dataSize(nrSession,2:3));
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
                        % calculate across-whole-session similarity values:
                        if wfr == 1
                            N = Priv.dataSize(nrSession,4);
                            %tic
                            if Pub.ssiOn && ~flag(2)
                                % Sign-similarity index:
                                dts = sign(diff(ts));
                                ssiData(yy,zz) = single(sum(abs(sum(dts,2)))/numel(dts));
                            end
                            if Pub.nmiOn && ~flag(3)
                                % Mutual information of binary data:
                                dts(dts == -1) = 0;
                                % call binary_mi.m:
                                nmiData(yy,zz) = single(binary_mi(dts));
                            end
                            if Pub.corOn && ~flag(1)
                                % correlation coefficient calculation (ref: Matlab corrcoef.m)
                                [n1,m1] = size(ts);
                                xc = ts - repmat(sum(ts)/n1,n1,1);  % Remove mean
                                c1 = (xc' * xc) / (n1-1); % calculate inner products
                                d1 = sqrt(diag(c1)); % sqrt first to avoid under/overflow
                                dd = d1*d1';
                                dd(1:m1+1:end) = diag(c1); % remove roundoff on diag
                                r1 = c1 ./ dd;
                                r1 = r1(INDS);
                                % save all correlations:
                                %corMatData(yy,zz,:) = r1;
                                
                                % save mean correlation:
                                corData(yy,zz) = single(mean(r1));
                                
                                if mod(iter,10000) == 0
                                    disp(['Intersubject correlation value ('...
                                        num2str(xx) ',' num2str(yy) ',' ...
                                        num2str(zz) '): ' num2str(corData(yy,zz))])
                                end
                            end
                            % Kendall's coefficient of concordance:
                            % tic
                            if Pub.kenOn && ~flag(4)
                                [~,R] = sort(ts);
                                RS = sum(R,2);
                                S = sum(RS.^2)-N*mean(RS).^2;
                                F = Priv.nrSubjects*...
                                    Priv.nrSubjects*(N*N*N-N);
                                kenData(yy,zz) = single(12*S/F);
                                %kenT = toc
                            end
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % calculate similarity values for time-frames:
                        
                        if Priv.nrTimeIntervals(nrSession) > 0
                            N = Pub.windowSize;
                            ts_win = ts(Priv.startInds{nrSession}...
                                (wfr):Priv.startInds{nrSession}...
                                (wfr)+Pub.windowSize-1,:);
                            if Pub.ssiOn && ~flag(2+4)
                                % Sign-similarity index:
                                dts = sign(diff(ts_win));
                                ssiDatawin(yy,zz,wfr) = single(sum(abs(sum(dts,2)))/numel(dts));
                            end
                            if Pub.nmiOn && ~flag(3+4)
                                % Mutual information of binary data:
                                dts(dts == -1) = 0;
                                nmiDatawin(yy,zz,wfr) = single(binary_mi(dts));
                            end
                            if Pub.corOn && ~flag(1+4)
                                % Mean of pairwise correlation (ref: Matlab corrcoef.m):
                                [n1,m1] = size(ts_win);
                                xc = ts_win - repmat(sum(ts_win)/n1,n1,1);  % Remove mean
                                c1 = (xc' * xc) / (n1-1); % calculate inner products
                                d1 = sqrt(diag(c1)); % sqrt first to avoid under/overflow
                                dd = d1*d1';
                                dd(1:m1+1:end) = diag(c1); % remove roundoff on diag
                                r1 = c1 ./ dd;
                                r1 = r1(INDS);
                                % save all correlations:
                                %corMatDatawin(yy,zz,:,wfr) = r1;
                                % save mean correlation:
                                corDatawin(yy,zz,wfr) = single(mean(r1));
                            end
                            if Pub.kenOn && ~flag(4+4)
                                % Kendall's coefficient of concordance:
                                [~,R] = sort(ts_win);
                                RS = sum(R,2);
                                S = sum(RS.^2)-N*mean(RS).^2;
                                F = Priv.nrSubjects*...
                                    Priv.nrSubjects*(N*N*N-N);
                                kenDatawin(yy,zz,wfr) = single(12*S/F);
                            end
                        end
                    end
                end
            end
        end
        
        % save results:
        %    tic % result saving time may vary considerably
        if Pub.ssiOn && ~flag(2)
            mMapResult.whole.([Priv.prefixFreqBand...
                num2str(nrBand)]).([Priv.prefixSession num2str(nrSession)]...
                ).ssi.Data.xyz(xx,:,:) = ssiData;
        end
        if Pub.nmiOn && ~flag(3)
            mMapResult.whole.([Priv.prefixFreqBand...
                num2str(nrBand)]).([Priv.prefixSession...
                num2str(nrSession)]).nmi.Data.xyz(xx,:,:) = nmiData;
        end
        if Pub.corOn && ~flag(1)
            mMapResult.whole.([Priv.prefixFreqBand...
                num2str(nrBand)]).([Priv.prefixSession...
                num2str(nrSession)]).cor.Data.xyz(xx,:,:) = corData;
        end
        if Pub.kenOn && ~flag(4)
            mMapResult.whole.([Priv.prefixFreqBand...
                num2str(nrBand)]).([Priv.prefixSession ...
                num2str(nrSession)]).ken.Data.xyz(xx,:,:)= kenData;
        end
        if nrFrames > 1
            for wfr = 1:size(ssiDatawin,3)
                if Pub.ssiOn && ~flag(2+4)
                    mMapResult.win.([Priv.prefixFreqBand...
                        num2str(nrBand)]).([Priv.prefixSession...
                        num2str(nrSession)]).ssi.Data(wfr).xyz(xx,:,:) = ssiDatawin(:,:,wfr);
                end
                if Pub.nmiOn && ~flag(3+4)
                    mMapResult.win.([Priv.prefixFreqBand...
                        num2str(nrBand)]).([Priv.prefixSession...
                        num2str(nrSession)]).nmi.Data(wfr).xyz(xx,:,:) = nmiDatawin(:,:,wfr);
                end
                if Pub.corOn && ~flag(1+4)
                    mMapResult.win.([Priv.prefixFreqBand...
                        num2str(nrBand)]).([Priv.prefixSession...
                        num2str(nrSession)]).cor.Data(wfr).xyz(xx,:,:) = corDatawin(:,:,wfr);
                end
                if Pub.kenOn && ~flag(4+4)
                    mMapResult.win.([Priv.prefixFreqBand...
                        num2str(nrBand)]).([Priv.prefixSession...
                        num2str(nrSession)]).ken.Data(wfr).xyz(xx,:,:) = kenDatawin(:,:,wfr);
                end
            end
        end
        %      toc
    end
end

lock_name=['calcSimMaps_',num2str(nrBand),'_',num2str(nrSession)];
if(freeToWrite('check',Pub.dataDestination,lock_name))
    
    % set maps non-writable:
    if Pub.corOn
        if ~flag(1)
            load([Pub.dataDestination 'memMaps'])
            memMaps.resultMap.whole.([Priv.prefixFreqBand...
                num2str(nrBand)]).([Priv.prefixSession...
                num2str(nrSession)]).cor.Writable = false;
            save([Pub.dataDestination 'memMaps'],'memMaps')
            clear memMaps
        end
        if Priv.nrTimeIntervals(nrSession) > 0 && ~flag(1+4)
            load([Pub.dataDestination 'memMaps'])
            memMaps.resultMap.win.([Priv.prefixFreqBand...
                num2str(nrBand)]).([Priv.prefixSession...
                num2str(nrSession)]).cor.Writable = false;
            save([Pub.dataDestination 'memMaps'],'memMaps')
            clear memMaps
        end
    end
    if Pub.kenOn
        if ~flag(4)
            load([Pub.dataDestination 'memMaps'])
            memMaps.resultMap.whole.([Priv.prefixFreqBand...
                num2str(nrBand)]).([Priv.prefixSession...
                num2str(nrSession)]).ken.Writable = false;
            save([Pub.dataDestination 'memMaps'],'memMaps')
            clear memMaps
        end
        if Priv.nrTimeIntervals(nrSession) > 0 && ~flag(4+4)
            load([Pub.dataDestination 'memMaps'])
            memMaps.resultMap.win.([Priv.prefixFreqBand...
                num2str(nrBand)]).([Priv.prefixSession...
                num2str(nrSession)]).ken.Writable = false;
            save([Pub.dataDestination 'memMaps'],'memMaps')
            clear memMaps
        end
    end
    if Pub.ssiOn
        if ~flag(2)
            load([Pub.dataDestination 'memMaps'])
            memMaps.resultMap.whole.([Priv.prefixFreqBand...
                num2str(nrBand)]).([Priv.prefixSession...
                num2str(nrSession)]).ssi.Writable = false;
            save([Pub.dataDestination 'memMaps'],'memMaps')
            clear memMaps
        end
        if Priv.nrTimeIntervals(nrSession) > 0 ~flag(4+2)
            
            load([Pub.dataDestination 'memMaps'])
            memMaps.resultMap.win.([Priv.prefixFreqBand...
                num2str(nrBand)]).([Priv.prefixSession...
                num2str(nrSession)]).ssi.Writable = false;
            save([Pub.dataDestination 'memMaps'],'memMaps')
            clear memMaps
        end
    end
    if Pub.nmiOn
        if ~flag(3)
            load([Pub.dataDestination 'memMaps'])
            memMaps.resultMap.whole.([Priv.prefixFreqBand...
                num2str(nrBand)]).([Priv.prefixSession...
                num2str(nrSession)]).nmi.Writable = false;
            save([Pub.dataDestination 'memMaps'],'memMaps')
            clear memMaps
        end
        if Priv.nrTimeIntervals(nrSession) > 0 && ~flag(3+4)
            load([Pub.dataDestination 'memMaps'])
            memMaps.resultMap.win.([Priv.prefixFreqBand...
                num2str(nrBand)]).([Priv.prefixSession...
                num2str(nrSession)]).nmi.Writable = false;
            save([Pub.dataDestination 'memMaps'],'memMaps')
            clear memMaps
        end
    end
    
    [~]=freeToWrite('release',Pub.dataDestination,lock_name);
end
showTime(0);
