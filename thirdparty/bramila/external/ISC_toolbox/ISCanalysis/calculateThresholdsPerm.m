function calculateThresholdsPerm(Params,nrSession,winOn)


combineNullVals(Params,nrSession,winOn);
generateLUT(Params,nrSession,winOn);

load([Params.PublicParams.dataDestination 'memMaps'])
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

for nrBand = 0:Params.PrivateParams.maxScale + 1
  disp(['Band: ' num2str(nrBand)])
  vals = [];
  if winOn
    for k = 1:Params.PrivateParams.nrTimeIntervals(nrSession)
      disp(['time interval: ' num2str(k)])
      I = memMaps.(Params.PrivateParams.resultMapName).win.([...
      Params.PrivateParams.prefixFreqBand num2str(nrBand)]).([Params.PrivateParams.prefixSession...
      num2str(nrSession)]).cor.Data(k).xyz(:,:,:);
      vals = [vals ; I(bmask)];
    end
  else
    I = memMaps.(Params.PrivateParams.resultMapName).whole.([Params.PrivateParams.prefixFreqBand...
    num2str(nrBand)]).([Params.PrivateParams.prefixSession...
    num2str(nrSession)]).cor.Data.xyz(:,:,:);
    vals = I(bmask);
  end
  evaluatePvals(Params,vals,nrSession,winOn,nrBand);
end
