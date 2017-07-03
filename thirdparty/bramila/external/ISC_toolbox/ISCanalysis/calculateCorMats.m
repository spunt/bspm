function calculateCorMats(Params,nrBand,nrSession)

% This function calculates and saves all pairwise intersubject correlation maps.                          %
% If the number of subjects is N, the size of the voxelwise correlation matrix
% becomes N*N. Only non-redundant correlations are stored, i.e., N(N-1)/2 correlations
% per voxel. Correlation values can be accessed through memMaps.cormatMap -field.
%
% Inputs:
% Params - struct containing all necessary parameters
% nrBand - frequency subband index (note: 0 refers to full band)
% nrSession - session index
%
% Example:
% Calculate pairwise correlations for subband 4, session 1:
% load('analysisParameters','Params')
% calculateCorMats(Params,4,1);
% Then access correlation coefficients of the first 10 time-intervals from one voxel (x=25,y=40,z=50):
% N = Params.PrivateParams.nrSubjects;CorVals = zeros(10,N*(N-1)/2);
% load memMaps
% for k = 1:10
%    CorVals(k,:) = squeeze(memMaps.cormatMap.win.band4.Session1.cor.(['timeInt' num2str(k)]).Data.xyzc(25,40,50,:));
% end
%
% See also:
% ISCANALYSIS
% RUNANALYSIS
% CALCULATESIMILARITYMAPS
%

% Last modified 5.8.2013 by Juha Pajula
% Tampere University of Technology
% Department of Signal Processing
% e-mail: juha.pajula@tut.fi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

showTime(1);

Pub = Params.PublicParams;
Priv = Params.PrivateParams;

if ~Pub.calcCorMatrices
    disp('Correlation matrix calculation not specified, calculation skipped...')
    return
end

if Pub.corOn == 0
    disp('Correlation similarity must be selected when calculating correlation matrices -> skipped...')
    return
end

% load memory maps:
load([Pub.dataDestination 'memMaps'])
if nrBand == 0
    mMapOrig = memMaps.(Priv.origMapName).([Priv.prefixSession num2str(nrSession)]);
else
    mMapFilt = memMaps.(Priv.filtMapName).([Priv.prefixSession num2str(nrSession)]);
end
mMapCorMat = memMaps.(Priv.cormatMapName);
clear memMaps

if mMapCorMat.whole.([Priv.prefixFreqBand...
        num2str(nrBand)]).([Priv.prefixSession...
        num2str(nrSession)]).cor.Writable == false
    disp('Correlation matrices written already, canceling computation...')
    return
end
for wfr = 1:Priv.nrTimeIntervals(nrSession)
    if mMapCorMat.win.([Priv.prefixFreqBand...
            num2str(nrBand)]).([Priv.prefixSession...
            num2str(nrSession)]).cor.([Priv.prefixTimeVal ...
            num2str(wfr)]).Writable == false
        disp('Correlation matrices written already, canceling computation...')
        return
    end
end






[~,~,en] = computer;

% init time-series data matrix (across-subject):
cDat = zeros([Priv.dataSize(nrSession,[4 2 3]), Priv.nrSubjects]);
% load brain Mask:
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

% INDS = [];
% for hh = 1:Priv.nrSubjects-1
%     INDS = [INDS (1+hh*Priv.nrSubjects):(...
%         1+hh*Priv.nrSubjects+hh-1)];
% end

INDS = find(triu(ones(Priv.nrSubjects,Priv.nrSubjects),1));

iter = 0;

if Priv.nrTimeIntervals(nrSession) == 0
    nrFrames = 1;
else
    nrFrames = Priv.nrTimeIntervals(nrSession);
end


%%%%%%%%%%%%%%%%%%%%
%warning('Modified window steps!!!!')
%clear sI;for k=1:24;if mod(k,5)==0;sI(k)=8;else;sI(k)=9;end;end;sI=cumsum(sI)-8
%Priv.startInds{1} = sI;
%clear sI;for k=1:39;if mod(k,5)==0;sI(k)=8;else;sI(k)=9;end;end;sI=cumsum(sI)-8
%Priv.startInds{1} = sI;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation of similarity maps:

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
        corMatData = zeros([Priv.dataSize(nrSession,2:3),length(INDS)]);
        corMatDatawin = zeros([Priv.dataSize(nrSession,2:3),length(INDS),Priv.nrTimeIntervals(nrSession)]);
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
                                % save all correlations:
                                corMatData(yy,zz,:) = r1;
                                if mod(iter,10000) == 0
                                    disp(['Intersubject correlation values ('...
                                        num2str(xx) ',' num2str(yy) ',' ...
                                        num2str(zz) '): ' num2str(corMatData(yy,zz,:))])
                                end
                            end
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % calculate similarity values for time-frames:
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
                                % save all correlations:
                                corMatDatawin(yy,zz,:,wfr) = r1;
                            end
                        end
                    end
                end
            end
        end
        
        % save results:
        %      tic
        if Pub.corOn
            mMapCorMat.whole.([Priv.prefixFreqBand...
                num2str(nrBand)]).([Priv.prefixSession...
                num2str(nrSession)]).cor.Data.xyzc(xx,:,:,:) = corMatData;
        end
        for wfr = 1:Priv.nrTimeIntervals(nrSession)
            if Pub.corOn
                mMapCorMat.win.([Priv.prefixFreqBand...
                    num2str(nrBand)]).([Priv.prefixSession...
                    num2str(nrSession)]).cor.([Priv.prefixTimeVal ...
                    num2str(wfr)]).Data.xyzc(xx,:,:,:) = corMatDatawin(:,:,:,wfr);
            end
        end
        %      toc
    end
end

lock_name=['calcCorMats_',num2str(nrBand),'_',num2str(nrSession)];
if(freeToWrite('check',Pub.dataDestination,lock_name))
    % set maps non-writable:
    load([Pub.dataDestination 'memMaps'])
    memMaps.(Priv.cormatMapName).whole.([Priv.prefixFreqBand...
        num2str(nrBand)]).([Priv.prefixSession...
        num2str(nrSession)]).cor.Writable = false;
    for wfr = 1:Priv.nrTimeIntervals(nrSession)
        memMaps.(Priv.cormatMapName).win.([Priv.prefixFreqBand...
            num2str(nrBand)]).([Priv.prefixSession...
            num2str(nrSession)]).cor.([Priv.prefixTimeVal ...
            num2str(wfr)]).Writable = false;
    end
    
    save([Pub.dataDestination 'memMaps'],'memMaps')
    clear memMaps
    
    [~]=freeToWrite('release',Pub.dataDestination,lock_name);
end

showTime(0);
