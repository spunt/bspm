function [PrivateParams,PublicParams] = setPrivParams(PublicParams)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRIVATE PARAMETERS:
% Do not change these parameters unless necessary.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PrivateParams.subjectDestination = [PublicParams.dataDestination 'fMRIpreprocessed/'];
PrivateParams.subjectFiltDestination = [PublicParams.dataDestination 'fMRIfiltered/'];
PrivateParams.resultsDestination = [PublicParams.dataDestination 'results/'];
PrivateParams.statsDestination = [PublicParams.dataDestination 'stats/'];
PrivateParams.PFDestination = [PublicParams.dataDestination 'PF/'];
PrivateParams.PFsessionDestination = [PublicParams.dataDestination 'PFsession/'];
PrivateParams.withinDestination = [PublicParams.dataDestination 'within/'];
PrivateParams.phaseDifDestination = [PublicParams.dataDestination 'phase/'];

%PrivateParams.logDestination = [PublicParams.dataDestination 'scripts/'];

% set prefixes for the data files:
PrivateParams.prefixResults = 'simMeasure';
PrivateParams.prefixSubject = 'fMRIpreproc';
PrivateParams.prefixSubjectFilt = 'fMRIfilt';
PrivateParams.prefixSyncResults = 'synch';
PrivateParams.prefixPhaseSyncResults = 'phaseSynch';
PrivateParams.prefixLUT = 'LUT';
PrivateParams.prefixTimeVal = 'timeInt';
PrivateParams.prefixPF = 'PF';
PrivateParams.prefixPFMat = 'PFMat';
PrivateParams.prefixSessComp = 'sessComp';
%PrivateParams.prefixPFMat = 'PFsessionMat';

PrivateParams.prefixCorMat = 'corMat';
PrivateParams.prefixFreqComp = 'freqComp';
PrivateParams.prefixTMap = 'tstats';
PrivateParams.prefixWithin = 'within';
PrivateParams.prefixPhaseDif = 'phase';

PrivateParams.prefixSession = 'Session';
PrivateParams.prefixFreqBand = 'band';

PrivateParams.simM = [{'ssi'},{'nmi'},{'cor'},{'ken'}];

PrivateParams.nrSubjects = size(PublicParams.subjectSource,2);
PrivateParams.nrSessions = size(PublicParams.subjectSource,1);
PrivateParams.transformType = 'DWT';
PrivateParams.resultMapName = 'resultMap';
PrivateParams.origMapName = 'origMap';
PrivateParams.filtMapName = 'filtMap';
PrivateParams.synchMapName = 'synchMap';
PrivateParams.phaseSynchMapName = 'phaseSynchMap';
PrivateParams.statMapName = 'statMap';
PrivateParams.cormatMapName = 'cormatMap';
PrivateParams.PFMapName = 'PFMap';
PrivateParams.PFmatMapName = 'PFmatMap';
PrivateParams.withinMapName = 'withinMap';
PrivateParams.phaseMapName = 'phaseMap';

PrivateParams.PFMapSessionName = 'PFMapSession';
PrivateParams.PFmatMapSessionName = 'PFMatMapSession';

if ~PublicParams.corOn
  PublicParams.calcStats = 0;
  PublicParams.calcCorMatrices = 0;
end

% get size of the data set:

% The following check is needed to avoid overwriting dataSize-field
% if it has been set in memMapdata.m. dataSize-field can be set in
% memMapdata.m if load_nii.m (called by memMapData) performs
% affine transformation, leading to mismatch between true data size
% and header information. Otherwise, header information is used to deter-
% mine the size (below).
if ~isfield(PrivateParams,'dataSizeMismatch')
    for m = 1:PrivateParams.nrSessions % all sessions
        fileFormat = PublicParams.fileFormatSubj;
        fileName = PublicParams.subjectSource{m,1};
        siz = getDataSize(fileName,fileFormat);
        PublicParams.dataSize(m,:) = siz;
    end
    PrivateParams.dataSize = PublicParams.dataSize;
end


for m = 1:PrivateParams.nrSessions % all sessions
  N = PrivateParams.dataSize(m,4);
  PrivateParams.nrTimeIntervals(m) = floor( (N-PublicParams.windowSize)...
  /PublicParams.windowStep ) + 1;
  PrivateParams.startInds{m} = 1:(PublicParams.windowStep):N;
  PrivateParams.startInds{m} = ...
  PrivateParams.startInds{m}(1:PrivateParams.nrTimeIntervals(m));
end


for m = 1:PrivateParams.nrSessions  
  if ~PublicParams.winOn
    PrivateParams.nrTimeIntervals(m) = 0;
  end
end


if ( PrivateParams.nrSessions ~= size(PrivateParams.dataSize,1) )
  error('Given data do not match the number of sessions!!')
  return
end

if PrivateParams.nrSessions > 1
    if PublicParams.sessionCompOn
        Len = PrivateParams.dataSize(:,end);
        if sum(diff(Len)) ~= 0
            error('To allow session comparison, time series length between sessions must be equal!')
            return
        end
    end
end
    
if PublicParams.useTemplate == 1
    PublicParams.maskPath = PublicParams.atlasPath;
    PublicParams.fileFormat = 'nii'; % standard template format
    if isequal(PrivateParams.dataSize(1,1:3),[182 218 182])
        PrivateParams.voxelSize = 1;
    elseif isequal(PrivateParams.dataSize(1,1:3),[91 109 91])
        PrivateParams.voxelSize = 2;
        PrivateParams.freqBlocks(1).block = 1;
        PrivateParams.freqBlocks(2).block = 2;
        PrivateParams.freqBlocks(3).block = 3;
        PrivateParams.freqBlocks(4).block = 4;
        PrivateParams.freqBlocks(5).block = 5;
        PrivateParams.freqBlocks(6).block = 6;
        PrivateParams.freqBlocks(7).block = 7;
    else
        error('Data size does not match known templates 91x109x91 (2mm) or 182x218x182 (1mm)!')
        return
    end

    PrivateParams.brainMask = [PublicParams.maskPath 'MNI152_T1_' ...
        num2str(PrivateParams.voxelSize) 'mm_brain_mask.nii'];
    atype = [{'cort'};{'sub'};{'cort'};{'sub'};{'cort'};{'sub'}];
    PrivateParams.atlasTh = [0 0 25 25 50 50];
    for k = 1:length(atype)
        PrivateParams.brainAtlases{k} = [PublicParams.atlasPath 'HarvardOxford-'...
            atype{k} '-maxprob-thr' num2str(PrivateParams.atlasTh(k)) '-' num2str(PrivateParams.voxelSize) 'mm.nii'];
    end

    % obtain brain regions:
    atD = load_nii([PrivateParams.brainAtlases{1}]);
    atD = atD.img;
    brainRegions = double([sort(nonzeros(unique(atD(:))))']);
    PrivateParams.nrRegions(1) = length(brainRegions);
    atD = load_nii([PrivateParams.brainAtlases{2}]);
    atD = atD.img;
    brainRegions2 = double([sort(nonzeros(unique(atD(:))))']);
    PrivateParams.nrRegions(2) = length(brainRegions2);
    PrivateParams.brainRegions = [brainRegions brainRegions2 255];
    clear atD;
else % custom binary mask
    PublicParams.maskPath = PublicParams.atlasPath;
    if length(PublicParams.maskPath) <= 4
        error('Full filename for binary mask must be given!')
    end    
    fileExt = PublicParams.maskPath(end-2:end);
    if strcmp(fileExt,'nii')
        Mask = load_nii(PublicParams.maskPath);
        Mask_dim=Mask.original.hdr.dime.dim;
        Mask_dim=Mask_dim(1:4);
        Mask = Mask.img;
        
    elseif strcmp(fileExt,'mat')
        Mask = load(PublicParams.maskPath);
        fiel = fields(Mask);
        Mask = Mask.(fiel{1});
        Mask = single(Mask);
        Mask_dim = [length(size(Mask)),size(Mask)];
    else
        error('Binary mask extension must be .nii or .mat!')
    end
    PublicParams.fileFormat = fileExt;
    if ~isnumeric(Mask)
        error('Binary mask must be numeric!')
    end
%    if length(size(Mask)) ~= 3
    if Mask_dim(1) ~= 3
        error('Binary mask must be 3-dimensional!')
    end
    if sum( Mask_dim(2:4) == PublicParams.dataSize(1,1:3) ) ~= 3
        error('Binary mask size does not fit fMRI data!')
    end
    VALS = unique(Mask(:));
    if length(VALS) ~= 2
        error('Mask is not binary!')
    end
    if sum( VALS(:) == [0 1]' ) ~= 2
        error('Binary mask must contain only ones and zeros!')
    end    
    PrivateParams.brainMask = PublicParams.maskPath;
end

%handles = setAtlasList(handles);


for k = 1:PrivateParams.nrSessions
  if PrivateParams.nrTimeIntervals(k) == 0 
    PrivateParams.winOn(k) = 0;
  else
    PrivateParams.winOn(k) = 1;
  end
end


% Set thresholds for synchronization curves:
PrivateParams.th =[0.05:0.05:0.5];

if PublicParams.nrFreqBands == 1
  PublicParams.nrFreqBands = 0;
end

% set maximum scale for SWT:
PrivateParams.maxScale = PublicParams.nrFreqBands - 1;

% control for modification of an existing project:
PrivateParams.filterUpdate = true;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Permutation test related fields:

% Set permutation test parameters. To save time and memory, calculation can be done parallel using 
% several CPU:s (nrPermutationSets). Total size of the null distribution used in statistical testing
% will be nrPermutationSets x nrPermutations.
PrivateParams.nrPermutationSets = PublicParams.nrPermutationSets;
PrivateParams.nrPermutations = PublicParams.nrPermutations;

% false-discovery rate before multiple corrections:
PrivateParams.q = 0.05;
% set null distribution look-up table intervals:
%PrivateParams.intVals = [0 0.005 0.01 0.02 0.03 0.04 0.05];
PrivateParams.intVals = [0 0.05];
% set approximate number of samples in each look-up table:
for k = 2:length(PrivateParams.intVals);
  PrivateParams.acc(k-1) = round( (PrivateParams.intVals(k)-PrivateParams.intVals(k-1)...
  )*PrivateParams.nrPermutationSets*PrivateParams.nrPermutations/(k-1));
end

if min(PrivateParams.acc) < 5
    error('Too few samples to conduct appropriate permutation tests!')
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computer information:

% Return information about computer on which MATLAB is running:
[str,maxsize,endian] = computer;
PrivateParams.computerInfo.endian = endian;
PrivateParams.computerInfo.str = 'str';
PrivateParams.computerInfo.maxsize = maxsize;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set paths:

% create data paths:
P.PrivateParams = PrivateParams;
P.PublicParams = PublicParams;

setDataPaths(P);

