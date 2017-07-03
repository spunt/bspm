function calculateSynchCurves(Params,nrBand,nrSession)

% This function calculates ROI-specific temporal intersubject synchronization (ISS) curves.
% Calcuted curves represent:
% 1) mean ISS in specific ROI
% 2) median ISS in specific ROI
% 3) number of voxels exceeding the threshold in specific ROI
%
% ISS calculations are performed across all frequency subbands.
% ROIs are based on Harvard-Oxford cortical and subcortical
% atlases (part of the FSL software package).
% Curves are calculated with 3 different atlas thresholds: 0%, 25%, and 50%.

% Last modified 5.8.2013 by Juha Pajula
% Tampere University of Technology 
% Department of Signal Processing
% e-mail: juha.pajula@tut.fi


showTime(1);

Pub = Params.PublicParams;
Priv = Params.PrivateParams;

if ~Pub.calcStandard
    disp('Basic ISC analysis not selected...')
    return
end

if ~Pub.useTemplate
    disp('Standard templates not in use, cannot compute ROI-based curves...')
    return
end

if Pub.calcPhase == 0 && Priv.nrTimeIntervals(nrSession) == 0
    disp('Time window ISC or ISPS not selected...')
    return
end

if Priv.nrTimeIntervals(nrSession) == 0
    corSynch = 0;
else
    corSynch = 1;
end

load([Pub.dataDestination 'memMaps.mat'])
[~,~,en] = computer;
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

flagSync = [false false];
if Pub.calcPhase
    if memMaps.(Priv.phaseSynchMapName).([Priv.prefixSession num2str(nrSession)]).([...
            Priv.prefixFreqBand num2str(nrBand)]).Writable == false
        disp('Phase synch curves for ROIs already computed...')
        flagSync(1) = true;
    end
else
    flagSync(1) = true;        
end
if Pub.winOn
    if memMaps.(Priv.synchMapName).([Priv.prefixSession num2str(nrSession)]).([...
            Priv.prefixFreqBand num2str(nrBand)]).Writable == false
        disp('Time window ISC curves for ROIs already computed...')
       flagSync(2) = true; 
    end
  
else
    flagSync(2) = true;
end
if isequal(flagSync,[true true])
        disp('Temporal ISC curves for ROIs already computed...')
    return
end

if corSynch
    load([Priv.statsDestination 'Thband' num2str(nrBand)...
            'Session' num2str(nrSession) 'win1'])
    LL = length(Th);
    % init synch.curve data matrix:
    T = zeros([Priv.nrTimeIntervals(nrSession),LL+2,...
        length(Priv.simM),length(Priv.brainAtlases)/2,length(Priv.brainRegions)]);

    % calculate synch.curves:
    atlasIter = 0;
    for u = 1:2:length(Priv.brainAtlases)
        disp(['Atlas: ' num2str(atlasIter) '/' ...
            num2str(length(Priv.brainAtlases)/2)])
        % obtain atlas:
        atCort = load_nii([Priv.brainAtlases{u}]);
        atlas = atCort.img;clear atCort;
        atSub = load_nii([Priv.brainAtlases{u+1}]);
        atlas(:,:,:,2) = atSub.img;clear atSub;
        atlasIter = atlasIter + 1;
        for n = 1:length(Priv.simM)
            if ( ( Pub.ssiOn && strcmp(Priv.simM{n},'ssi') ) || ...
                    ( Pub.nmiOn && strcmp(Priv.simM{n},'nmi') ) || ...
                    ( Pub.corOn && strcmp(Priv.simM{n},'cor') ) || ...
                    ( Pub.kenOn && strcmp(Priv.simM{n},'ken') ) )
                disp(['  Similarity measure: ' num2str(n) '/' ...
                    num2str(length(Priv.simM))])
                for t = 1:Priv.nrTimeIntervals(nrSession)
                    disp(['    Time interval: ' num2str(t) '/' ...
                        num2str(Priv.nrTimeIntervals(nrSession))])
                    if ~strcmp(Priv.computerInfo.endian,en)
                        R = swapbytes(memMaps.(Priv.resultMapName).win.([Priv.prefixFreqBand...
                            num2str(nrBand)]).([Priv.prefixSession num2str(nrSession)]).(Priv.simM{n}).Data(t).xyz);
                    else
                        R = memMaps.(Priv.resultMapName).win.([Priv.prefixFreqBand...
                            num2str(nrBand)]).([Priv.prefixSession num2str(nrSession)]).(Priv.simM{n}).Data(t).xyz;
                    end
                    for p = 1:length(Priv.brainRegions)
                        % voxel amount based curves:
                        for s = 1:LL
                            if Priv.brainRegions(p) == 255 % whole brain synch. curves
                                T(t,s,n,atlasIter,p) = sum(R(find(bmask & ~isnan(R))) >= Th(s));
                            else % ROI synch. curves
                                if p <= Priv.nrRegions(1)
                                    atInd = 1;
                                else
                                    atInd = 2;
                                end
                                T(t,s,n,atlasIter,p) = sum( R(find( atlas(:,:,:,atInd) == ...
                                    Priv.brainRegions(p) & ~isnan(R) ) ) >= Th(s));
                            end
                        end
                        % mean and median curves
                        if Priv.brainRegions(p) == 255 % whole brain
                            T(t,LL+1,n,atlasIter,p) = mean(R(find(bmask & ~isnan(R))));
                            T(t,LL+2,n,atlasIter,p) = median(R(find(bmask & ~isnan(R))));
                        else % ROI curves:
                            if p <= Priv.nrRegions(1)
                                atInd = 1;
                            else
                                atInd = 2;
                            end
                            T(t,LL+1,n,atlasIter,p) = mean(R(find( atlas(:,:,:,atInd) == ...
                                Priv.brainRegions(p) & ~isnan(R) )));
                            T(t,LL+2,n,atlasIter,p) = median(R(find( atlas(:,:,:,atInd) == ...
                                Priv.brainRegions(p) & ~isnan(R) )));
                        end
                    end
                end
            end
        end
    end
    
    if memMaps.(Priv.synchMapName).([Priv.prefixSession num2str(nrSession)]).([...
            Priv.prefixFreqBand num2str(nrBand)]).Writable == true
        
        lock_name=['calcSynchCurv_',num2str(nrBand),'_',num2str(nrSession)];
        if(freeToWrite('check',Pub.dataDestination,lock_name))
            for p = 1:length(Priv.brainRegions)
                memMaps.(Priv.synchMapName).([Priv.prefixSession num2str(nrSession)]).([...
                    Priv.prefixFreqBand num2str(nrBand)]).Data(p).tcsa = single(T(:,:,:,:,p));
            end
            memMaps.(Priv.synchMapName).([Priv.prefixSession num2str(nrSession)]).([...
                Priv.prefixFreqBand num2str(nrBand)]).Writable = false;
            save([Pub.dataDestination 'memMaps'],'memMaps')
            [~]=freeToWrite('release',Pub.dataDestination,lock_name);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Pub.calcPhase == 1

    % init synch.curve data matrix:
    Tp = zeros([Priv.dataSize(nrSession,4),2,length(Priv.brainAtlases)/2,length(Priv.brainRegions)]);
    B = construct4D(Params,nrSession,nrBand,2);
    % calculate synch.curves:
    atlasIter = 0;
    for u = 1:2:length(Priv.brainAtlases)
        disp(['Atlas Type: ' num2str(atlasIter) '/' ...
            num2str(length(Priv.brainAtlases)/2)])
        % obtain atlas:
        atCort = load_nii([Priv.brainAtlases{u}]);
        atlas = atCort.img;clear atCort;
        atSub = load_nii([Priv.brainAtlases{u+1}]);
        atlas(:,:,:,2) = atSub.img;clear atSub;
        atlasIter = atlasIter + 1;

        for t = 1:Priv.dataSize(nrSession,4)
            Bt = B(:,:,:,t);
            for p = 1:length(Priv.brainRegions)
                % mean and median curves
                if Priv.brainRegions(p) == 255 % whole brain
                    Tp(t,1,atlasIter,p) = mean(Bt(find(bmask & ~isnan(Bt))));
                    Tp(t,2,atlasIter,p) = median(Bt(find(bmask & ~isnan(Bt))));
                else % ROI curves:
                    if p <= Priv.nrRegions(1)
                        atInd = 1;
                    else
                        atInd = 2;
                    end
                    Tp(t,1,atlasIter,p) = mean(Bt(find( atlas(:,:,:,atInd) == ...
                        Priv.brainRegions(p) & ~isnan(Bt) )));
                    Tp(t,2,atlasIter,p) = median(Bt(find( atlas(:,:,:,atInd) == ...
                        Priv.brainRegions(p) & ~isnan(Bt) )));
                end
            end
        end
    end
    if memMaps.(Priv.phaseSynchMapName).([Priv.prefixSession num2str(nrSession)]).([...
            Priv.prefixFreqBand num2str(nrBand)]).Writable == true
        
        lock_name=['calcSynchCurv_',num2str(nrBand),'_',num2str(nrSession)];
        if(freeToWrite('check',Pub.dataDestination,lock_name))
            for p = 1:length(Priv.brainRegions)
                memMaps.(Priv.phaseSynchMapName).([Priv.prefixSession num2str(nrSession)]).([...
                    Priv.prefixFreqBand num2str(nrBand)]).Data(p).tca(:,1:2,:) = single(Tp(:,:,:,p));
            end
            memMaps.(Priv.phaseSynchMapName).([Priv.prefixSession num2str(nrSession)]).([...
                Priv.prefixFreqBand num2str(nrBand)]).Writable = false;
            save([Pub.dataDestination 'memMaps'],'memMaps')
            [~]=freeToWrite('release',Pub.dataDestination,lock_name);
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



showTime(0);
