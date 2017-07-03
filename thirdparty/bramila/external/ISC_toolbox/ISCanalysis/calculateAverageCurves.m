function calculateAverageCurves(Params,nrSubject)

% calculate regional average time-series of a subject
% 
% 6.10.2009

Priv = Params.PrivateParams;
Pub = Params.PublicParams;

disp(['Calculate regional average time-series for subject ' num2str(nrSubject)])

for nrSession = 1:Priv.nrSessions
  meanTS = zeros(Priv.maxScale+2,3,(length(Priv.brainRegions)-1),Priv.dataSize(nrSession,4));
  size(meanTS)
  disp(['session ' num2str(nrSession)])
  for nrBand = 0:Priv.maxScale + 1;
    disp([' band ' num2str(nrBand)])
    % construct 4D fMRI data from x-blocks:
    B = construct4D(Params,nrSession,nrBand,nrSubject);
    SIZ = size(B);
    atlasIter = 0;
    for u = 1:2:length(Priv.brainAtlases)
      atlasIter = atlasIter + 1;
      % obtain atlas:
      atCort = load_nii([Priv.brainAtlases{u}]);
      atlas = atCort.img;clear atCort;
      atSub = load_nii([Priv.brainAtlases{u+1}]);
      atlas(:,:,:,2) = atSub.img;clear atSub;
      Priv.brainAtlases{u}
      for atlasRegion = 1:(length(Priv.brainRegions)-1)
        % switch from cortical to subcortical atlas when all
        % cortical ROIs have been calculated:
        if atlasRegion <= Priv.nrRegions(1)
          atInd = 1;
        else
          atInd = 2;
        end
        disp(['     atlas region ' num2str(atlasRegion)])
        % find indices of a ROI:
        atInds = find(atlas(:,:,:,atInd) == Priv.brainRegions(atlasRegion));
        disp(['atlas region: ' num2str(Priv.brainRegions(atlasRegion))])
        % obtain all time-series of the ROI:
        TS = zeros(prod(SIZ(1:3)),size(B,4));
        for k = 1:size(B,4)
          X = B(:,:,:,k);
          TS(:,k) = X(:);
        end
        TS = TS(atInds,:);
        rowSum = sum(TS,2);
        % remove nans if any:
        nans = isnan(rowSum);
        sum(nans)
        TS(nans,:) = [];
        % calculate mean time-series:
        mTS = mean(TS);
        meanTS(nrBand+1,atlasIter,atlasRegion,:) = mTS;
      end
    end
  end
  save([Priv.withinDestination 'meanTimeSeriesSession' num2str(nrSession) 'Subject' num2str(nrSubject)],'meanTS')
end
