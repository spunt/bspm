function calculateWithinSubjectCorrelation(Params,nrBand,nrSession)

% This function calculates intersubject synchronization maps 
% according to user specified parameters.
% 
% Inputs:
% Params - analysis parameters from initParams.m
% nrBand - frequency subband index (0 = full, original band)
% nrSession - session index
%

% example:
% load('/destination/crash_testi.mat')
% calculateWithinSubjectCorrelation(Params,0,1);

%
% See also:
% INITPARAMS
% RUNanalysis

% 19.09.2009
% Jukka-Pekka Kauppi
% Tampere University of Technology


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pub = Params.PublicParams;
Priv = Params.PrivateParams;
load([Pub.dataDestination 'memMaps'])
dMap.(Priv.origMapName) = memMaps.(Priv.origMapName);
dMap.(Priv.filtMapName) = memMaps.(Priv.filtMapName);
mMap.(Priv.withinMapName) = memMaps.(Priv.withinMapName);
clear memMaps
[pl,ms,en] = computer;

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

tsROI = zeros(Priv.dataSize(nrSession,4),Priv.nrSubjects);
size(tsROI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation of similarity maps:
for ROI_ind = 1:(length(Priv.brainRegions)-1)
  for subj = 1:Priv.nrSubjects
    load([Priv.withinDestination 'meanTimeSeriesSession' num2str(nrSession) 'Subject' num2str(subj)])
    mTS = squeeze(meanTS(nrBand+1,3,ROI_ind,:));
    tsROI(:,subj) = mTS;
  end
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
            cDat(:,:,:,k) = swapbytes(dMap.(Priv.origMapName).([Priv.prefixSession ...
            num2str(nrSession)]).([Priv.prefixSubject ...
            num2str(k)]).Data(xx).tyz);
          else
            cDat(:,:,:,k) = dMap.(Priv.origMapName).([Priv.prefixSession ...
            num2str(nrSession)]).([Priv.prefixSubject ...
            num2str(k)]).Data(xx).tyz(:,:,:);
          end
        else % load sub-band data
          if ~strcmp(Priv.computerInfo.endian,en)
            cDat(:,:,:,k) = swapbytes(dMap.(Priv.filtMapName).([Priv.prefixSession ...
            num2str(nrSession)]).([Priv.prefixSubjectFilt ...
            num2str(k)]).([Priv.prefixFreqBand num2str(nrBand)]).Data(xx).tyz);
          else
            cDat(:,:,:,k) = dMap.(Priv.filtMapName).([Priv.prefixSession ...
            num2str(nrSession)]).([Priv.prefixSubjectFilt ...
            num2str(k)]).([Priv.prefixFreqBand num2str(nrBand)]).Data(xx).tyz;
          end
        end
      end      
      N = Priv.dataSize(nrSession,4);
      corData = zeros([Priv.dataSize(nrSession,2:3) Priv.nrSubjects]);
      for yy = 1:Priv.dataSize(nrSession,2)
        for zz = 1:Priv.dataSize(nrSession,3)     
          if bmask(xx,yy,zz)
            % obtain each subject's time series:
            ts = squeeze(cDat(:,yy,zz,:));
            for k = 1:Priv.nrSubjects
              TS = [ts(:,k) tsROI(:,k)];
              % correlation coefficient calculation (ref: Matlab corrcoef.m)
              [n1,m1] = size(TS);
              xc = TS - repmat(sum(TS)/n1,n1,1); % Remove mean
              c1 = (xc' * xc) / (n1-1); % calculate inner products
              d1 = sqrt(diag(c1)); % sqrt first to avoid under/overflow
              dd = d1*d1';
              dd(1:m1+1:end) = diag(c1); % remove roundoff on diag
              r1 = c1 ./ dd;
              r1 = r1(1,2);
              corData(yy,zz,k) = single(r1);
            end
          end
        end
      end
%      assignin('base','corDatayyzzk',corData)
%      assignin('base','TS_ts_roiTs',TS)
%      xx
%      ROI_ind
      % save results:
%      AA=mMap.(Priv.withinMapName).whole.([Priv.prefixFreqBand...
%      num2str(nrBand)]).([Priv.prefixSession...
%      num2str(nrSession)]).cor.Data(ROI_ind).xyzs(xx,:,:,:);
      mMap.(Priv.withinMapName).whole.([Priv.prefixFreqBand...                                                                                                                                             
      num2str(nrBand)]).([Priv.prefixSession...                                                                                                                                                               
      num2str(nrSession)]).cor.Data(ROI_ind).xyzs(xx,:,:,:) = corData;
%      BB=mMap.(Priv.withinMapName).whole.([Priv.prefixFreqBand...                                                                                                                                             
%      num2str(nrBand)]).([Priv.prefixSession...                                                                                                                                                               
%      num2str(nrSession)]).cor.Data(ROI_ind).xyzs(xx,:,:,:);
%      assignin('base','bef',AA)
%      assignin('base','aft',BB)
%    return
    end
  end
end
