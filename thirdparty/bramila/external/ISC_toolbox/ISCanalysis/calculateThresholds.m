function calculateThresholds(Params,nrSession,winOn,mapType)

% This function thresholds the mean inter-subject correlation maps (GISC maps)
% based on generated permutation distribution. The following thresholds are obtained:
%
% 1.  p < 0.05,  no multiple comparisons correction
% 2.  p < 0.05,  FDR corrected using independence or positive dependence assumption
% 3.  p < 0.05,  FDR corrected (no assumptions)
% 4.  p < 0.05,  Bonferroni corrected
% 5.  p < 0.005,  no multiple comparisons correction
% 6.  p < 0.005, FDR corrected using independence or positive dependence assumption
% 7.  p < 0.005, FDR corrected (no assumptions)
% 8.  p < 0.005, Bonferroni corrected
% 9.  p < 0.001,  no multiple comparisons correction
% 10. p < 0.001, FDR corrected using independence or positive dependence assumption
% 11. p < 0.001, FDR corrected (no assumptions)
% 12. p < 0.001, Bonferroni corrected
%
% inputs:
% Params - struct containing all necessary parameters
% nrSession - session index
% winOn - perform calculations for windowed (winOn=1) or across-session data (winOn=0)
%
% See also:
%
% COMBINENULLVALS
% GENERATELUT
% EVALUATEPVALS
% INITPARAMS
% RUNANALYSIS

% Last modified 06.10.2010 by Jukka-Pekka Kauppi
% Tampere University of Technology
% Department of Signal Processing
% e-mail: jukka-pekka.kauppi@tut.fi

if ~Params.PublicParams.calcStandard
    disp('Standard inter-subject synchronization analysis not specified, critical values not calculated....')
    return
end

if ~Params.PublicParams.winOn && winOn
    disp('Time-windows not used, skipping significance test...')
    return
end

% for permSetIdx = 1:Params.PrivateParams.nrPermutationSets
%   if winOn
%       Fi = [Params.PrivateParams.statsDestination 'vals' num2str(permSetIdx) 'win.mat'];
%     if exist(Fi) == 2
%         disp(['File: ' Fi ' already found, skipping significance test...'])
%         return
%     end
%   else
%       Fi = [Params.PrivateParams.statsDestination 'vals' num2str(permSetIdx) ...
%             Params.PrivateParams.prefixSession num2str(nrSession) '.mat'];
%     if exist(Fi) == 2
%         disp(['File: ' Fi ' already found, skipping significance test...'])
%         return
%     end
%   end
% end

[~,~,en] = computer;

combineNullVals(Params,nrSession,winOn); % combine (right-tail) null-distribution data
generateLUT(Params,nrSession,winOn); % create look-up tables

load([Params.PublicParams.dataDestination 'memMaps']) % load memory map objects
% load brain mask:
if strcmp(Params.PublicParams.fileFormat,'nii')
    bmask = load_nii(Params.PrivateParams.brainMask);
    bmask = single(bmask.img);
elseif strcmp(Params.PublicParams.fileFormat,'mat')
    bmask = load(Params.PrivateParams.brainMask);
    fiel = fields(bmask);
    bmask = bmask.(fiel{1});
    bmask = single(bmask);
else
    error('Mask must be mat- or nii-file!')
end
bmask = logical(bmask);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate thresholds for each frequency band:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% obtain data for which thresholds will be calculated:
if winOn
    mapType = 1;
else
    mapType = 1;
    if Params.PublicParams.calcStats
        mapType = [mapType 2];
    end
    if Params.PublicParams.calcPhase
        mapType = [mapType 3];
    end
end

% get observed data depending on the input "mapType":
for bb = mapType
    for nrBand = 0:Params.PrivateParams.maxScale + 1
%        disp(['Band: ' num2str(nrBand)])
        vals = [];
        if winOn % obtain windowed data:
            disp(['Get thresholds for time-window ISC: band ' num2str(nrBand)])
            for k = 1:Params.PrivateParams.nrTimeIntervals(nrSession)
                disp(['time interval: ' num2str(k)])
                I = memMaps.(Params.PrivateParams.resultMapName).win.([...
                    Params.PrivateParams.prefixFreqBand num2str(nrBand)]).([Params.PrivateParams.prefixSession...
                    num2str(nrSession)]).cor.Data(k).xyz(:,:,:);
                vals = [vals ; I(bmask)];
            end
        else % obtain across-session data:
            switch bb
                case 1
                    disp(['Get thresholds for ISC maps: band ' num2str(nrBand)])
                    I = memMaps.(Params.PrivateParams.resultMapName).whole.([Params.PrivateParams.prefixFreqBand...
                        num2str(nrBand)]).([Params.PrivateParams.prefixSession...
                        num2str(nrSession)]).cor.Data.xyz(:,:,:);
                    vals = I(bmask);
                case 2
                    disp(['Get thresholds for t-stat maps: band ' num2str(nrBand)])
                    I = memMaps.(Params.PrivateParams.statMapName).whole.([Params.PrivateParams.prefixFreqBand...
                        num2str(nrBand)]).([Params.PrivateParams.prefixSession...
                        num2str(nrSession)]).cor.Data.xyz(:,:,:);
                    vals = I(bmask);
                case 3
                    disp(['Get thresholds for ISPS maps: band ' num2str(nrBand)])
                    vals = [];
                    for jj = 1:Params.PrivateParams.dataSize(nrSession,1)
                        Im = memMaps.(Params.PrivateParams.phaseMapName).([Params.PrivateParams.prefixSession...
                            num2str(nrSession)]).([Params.PrivateParams.prefixFreqBand...
                            num2str(nrBand)]).Data(jj).tyz(:,:,:);
                        SS = squeeze(bmask(jj,:,:));
                        SSS = zeros(size(Im,1),size(SS,1),size(SS,2));
                        for mm = 1:size(SSS,1)
                            SSS(mm,:,:) = SS;
                        end
                        vals = [vals; Im(logical(SSS))];
                    end
            end
        end
        % swap bytes if needed:
        if ~strcmp(Params.PrivateParams.computerInfo.endian,en)
            I = swapbytes(I);
        end
        % calculate critical thresholds:
        evaluatePvals(Params,vals,nrSession,winOn,nrBand,bb);
        
    end
end


