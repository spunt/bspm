function TS = getTS(Params,nrSession,nrBand,atlasRegion,atlasType,atlasTh,winOn)
  % Obtain time-series showing statistically significant inter-subject
  % synchronization.
  %
  % inputs:
  % Params - LSA parameter-struct 
  % nrSession - session index
  % nrBand - frequency band index
  % atlasRegion - atlas region index
  % atlasType - type of atlas: 1=cortical, 2=subcortical
  % atlasTh - atlas threshold: 1=0%, 2=25%, 3=50%
  % winOn - load windowed data: 1=yes, 0=no
  %
  % output:
  % TS - time-series: TS{time-interval}(voxel,time-point,subject)
  %
  % EX.
  % TS = getTS(Params,1,5,1,1,3,1);


Priv = Params.PrivateParams;
Pub = Params.PublicParams;
S = 1:2:5;
if atlasType == 2
  S = S + 1;
end
Priv.brainAtlases{S(atlasTh)}
at = load_nii(Priv.brainAtlases{S(atlasTh)});
at = find(at.img == atlasRegion);
load([Pub.dataDestination 'memMaps'])
[pl,ms,en] = computer;

if winOn
  fn = 'R005';
else
  fn = 'r005';
end

R = load([Priv.resultsDestination 'Regions'],fn)

if winOn
  intVals = Priv.nrTimeIntervals(nrSession);
  TS = cell(Priv.nrTimeIntervals(nrSession),1);
else
  intVals = 1
  TS = cell(1,1);
end

for ti = 1:intVals
  disp(['Time interval: ' num2str(ti) '/' num2str(intVals)])
  if winOn
    linInds = R.(fn){nrSession,nrBand+1,ti};
    startInd = Priv.startInds{nrSession}(ti);
    stopInd = Priv.startInds{nrSession}(ti)+Pub.windowSize-1;
    N = Pub.windowSize;
  else
    linInds = R.(fn){nrSession,nrBand+1};
    startInd = 1;
    stopInd = Priv.dataSize(nrSession,4);
    N = stopInd;
  end
  linInds = intersect(linInds,at);
  [xx,yy,zz] = ind2sub(Priv.dataSize(nrSession,1:3),linInds);
  [xx sidx] = sort(xx);
  if isempty(xx)
    disp('No active voxels --> reading terminated')
    return
  end
  lastX = xx(end);
  yy = yy(sidx);
  zz = zz(sidx);
  linInds = linInds(sidx);
  
  A = zeros([N,Priv.dataSize(nrSession,2:3),Priv.nrSubjects]);
  ts = zeros(length(linInds),N,Priv.nrSubjects);
%  size(ts)
  iters = length(linInds)
  for iter = 1:iters
    disp(['x: ' num2str(xx(iter)) '/' num2str(lastX)])
    if iter == 1 || xx(iter) ~= xx(iter-1)
      for m = 1:Priv.nrSubjects
        if nrBand == 0
          if ~strcmp(en,Priv.computerInfo.endian)
            A(:,:,:,m) = swapbytes(memMaps.origMap.([Priv.prefixSession num2str(nrSession)]...
            ).([Priv.prefixSubject num2str(m)]).Data(xx(iter)).tyz(startInd:stopInd,:,:));%(:,yy(iter),zz(iter));
          else
            A(:,:,:,m) = memMaps.origMap.([Priv.prefixSession num2str(nrSession)]...
            ).([Priv.prefixSubject num2str(m)]).Data(xx(iter)).tyz(startInd:stopInd,:,:);%(:,yy(iter),zz(iter));
          end
        else
          if ~strcmp(en,Priv.computerInfo.endian)
            A(:,:,:,m) = swapbytes(memMaps.filtMap.([Priv.prefixSession num2str(nrSession)]...
            ).([Priv.prefixSubjectFilt num2str(m)]).([Priv.prefixFreqBand ...
            num2str(nrBand)]).Data(xx(iter)).tyz(startInd:stopInd,:,:));%(:,yy(iter),zz(iter));
          else
            A(:,:,:,m) = memMaps.filtMap.([Priv.prefixSession num2str(nrSession)]...
            ).([Priv.prefixSubjectFilt num2str(m)]).([Priv.prefixFreqBand ...
            num2str(nrBand)]).Data(xx(iter)).tyz(startInd:stopInd,:,:);%(:,yy(iter),zz(iter));
          end
          %          
          ts_win = ts(Priv.startInds{nrSession}(wfr):Priv.startInds{nrSession}(wfr)+Pub.windowSize-1,:);
        end
      end
    end
    ts(iter,:,:) = squeeze(A(:,yy(iter),zz(iter),:));
%    max(ts(:))
  end
  TS{ti} = ts;
  INDS{ti} = linInds;
end

save([Priv.tsDestination 'TimeSeriesDataSession' num2str(nrSession) ...
'Band' num2str(nrBand) 'Region' num2str(atlasRegion) 'atlasType' ...
num2str(atlasType) 'Th' num2str(atlasTh) 'win' num2str(winOn)],'TS','INDS')
