function calculatePhaseSynch(Params,nrBand,nrSession)

% This function calculates inter-subject phase-synchronization.
%
% Inputs:
% Params - struct containing all necessary parameters
% nrBand - frequency subband index (note: 0 refers to full band)
% nrSession - session index
%
% See also:
% RUNANALYSIS
% ISCANALYSIS
% INITPARAMS

% Last modified 20.11.2013 by Jukka-Pekka Kauppi
% University of Helsinki
% Department of Computer Science
% e-mail: jukka-pekka.kauppi@helsinki.fi
%
% Modified 5.8.2013 by Juha Pajula
% Tampere University of Technology
% Department of Signal Processing
% e-mail: juha.pajula@tut.fi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

showTime(1);

Pub = Params.PublicParams;
Priv = Params.PrivateParams;

if ~Pub.calcPhase
    disp('Phase inter-subject synchronization analysis not selected, calculation skipped...')
    return
end

% load mempory map objects:
load([Pub.dataDestination 'memMaps'])
if Pub.nrFreqBands > 0
    mMapFilt = memMaps.(Priv.filtMapName).([Priv.prefixSession num2str(nrSession)]);
end
mMapOrig = memMaps.(Priv.origMapName).([Priv.prefixSession num2str(nrSession)]);
mMapPhase = memMaps.(Priv.phaseMapName);
clear memMaps
try
    if ~isempty(mMapPhase)
        if mMapPhase.([Priv.prefixSession...
                num2str(nrSession)]).([Priv.prefixFreqBand...
                num2str(nrBand)]).Writable == false
            disp('Phase maps written already, calculation canceled...')
            return
        end
    else
        disp('Memory pointers not available, calculation canceled...')
        return
    end
catch err
    error('Cannot update phase map data possibly due to changes in analysis parameters. To update the data, delete bin-files in phase-folder and try analysis again!')
end

% get platform information:
[~,~,en] = computer;

N = Priv.dataSize(nrSession,4);

cDat = zeros([Priv.dataSize(nrSession,[4 2 3]), Priv.nrSubjects]);
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
    return
end
bmask = logical(bmask);

INDS = find(triu(ones(Priv.nrSubjects,Priv.nrSubjects),1));
subjPairs = length(INDS);

iter = 0;
tic

% for each voxel inside the brain:
for xx = 1:Priv.dataSize(nrSession,1)
    disp(['x:' num2str(xx) '/' num2str(Priv.dataSize(nrSession,1))])
    % process only non-zero slices:
    if sum(sum(squeeze(bmask(xx,:,:)))) > 0
        cDat = zeros(Priv.dataSize(nrSession,4),Priv.dataSize(nrSession,2),...
            Priv.dataSize(nrSession,3),Priv.nrSubjects);
        
        % get source time-series of the subjects:
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
        
        % init temp data matrix:
        ph = zeros(Priv.dataSize(nrSession,[4 2 3]));
        phMax = ph;
        
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
                    ts = squeeze(cDat(:,yy,zz,:))';
                    ts = ts - mean(ts,2)*ones(size(ts,2),1)';
                    % size(ts)                    
                    
                    Rm = zeros(subjPairs,N);
                    ite = 1;
                    for kk = 1:Priv.nrSubjects
                        y1 = hilbert(ts(kk,:));
                        for hh = 1:Priv.nrSubjects
                            if hh > kk
                                y2 = hilbert(ts(hh,:));
                                Rm(ite,:) = angle(y1.*conj(y2));
                                ite = ite + 1;
                            end
                        end
                    end
                    if subjPairs == 1
                        ph(:,yy,zz) = abs(Rm);
                        ph(:,yy,zz) = 1 - ph(:,yy,zz)/pi;
                        phMax(:,yy,zz) = ph(:,yy,zz);
                    else
                        ph(:,yy,zz) = sum(abs(Rm))/subjPairs;
                        ph(:,yy,zz) = 1 - ph(:,yy,zz)/pi;
                        phMax(:,yy,zz) = min(abs(Rm));
                        phMax(:,yy,zz) = 1 - phMax(:,yy,zz)/pi;
                    end                    
                    
                    if mod(iter,10000) == 0
                        PP = ph(:,yy,zz);
                        PP = PP(:)';
                        PP2 = phMax(:,yy,zz);
                        PP2 = PP2(:)';
                        disp(['phase difference values ('...
                            num2str(xx) ',' num2str(yy) ',' ...
                            num2str(zz) '): ' num2str(PP)])
                        disp(num2str(PP2))
                    end
                end
            end
        end
        %    size(ph)
        % write data to disk:
        mMapPhase.([Priv.prefixSession...
            num2str(nrSession)]).([Priv.prefixFreqBand...
            num2str(nrBand)]).Data(xx).tyz(:,:,:) = ph;
        
    end
end

lock_name=['calcPhaseSynch_',num2str(nrBand),'_',num2str(nrSession)];
if(freeToWrite('check',Pub.dataDestination,lock_name))
    load([Pub.dataDestination 'memMaps'])
    memMaps.(Priv.phaseMapName).([Priv.prefixSession...
        num2str(nrSession)]).([Priv.prefixFreqBand...
        num2str(nrBand)]).Writable = false;
    save([Pub.dataDestination 'memMaps.mat'],'memMaps')
    clear memMaps mMapPhase
    [~]=freeToWrite('release',Pub.dataDestination,lock_name);
end
