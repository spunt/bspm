function returnTimeSeries(Params,nrSession,regionIndex,atlasType,atlasTh)

% calculate regional average time-series of a subject
% 
% 6.10.2009

Priv = Params.PrivateParams;
Pub = Params.PublicParams;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select brain atlas:
S = 1:2:5;
if atlasType == 2
  S = S + 1;
end
Priv.brainAtlases{S(atlasTh)}
at = load_nii(Priv.brainAtlases{S(atlasTh)});
atInds = find(at.img == regionIndex);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['atlas region ' num2str(regionIndex)])
size(atInds)
for nrSubject = 1:Priv.nrSubjects
  disp(['session ' num2str(nrSession)])
  for nrBand = 0:Priv.maxScale + 1;
    disp([' band ' num2str(nrBand)])
    % construct 4D fMRI data from x-blocks:
    B = construct4D(Params,nrSession,nrBand,nrSubject);
    SIZ = size(B);
    % obtain all time-series of the ROI:
    TS = single( zeros(prod(SIZ(1:3)),size(B,4)) );
    for k = 1:size(B,4)
      X = B(:,:,:,k);
      TS(:,k) = single( X(:) );
    end
    TS = single( TS(atInds,:) );
    save([Pub.dataDestination 'tsSession' num2str(nrSession) 'Subject' ...
    num2str(nrSubject) 'Band' num2str(nrBand) 'atlasType' num2str(atlasType) ...
    'AtlasTh' num2str(atlasTh) 'Region' num2str(regionIndex)],'TS')
  end
end
