function B = construct4D(Params,nrSession,nrBand)
  % Obtain time-series showing statistically significant inter-subject
  % synchronization.
  %
  % inputs:
  % Params - parameter-struct 
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

load([Pub.dataDestination 'memMaps'])
[pl,ms,en] = computer;

if winOn
  intVals = Priv.nrTimeIntervals(nrSession);
else
  intVals = 1
end

B = zeros([Priv.dataSize(nrSession,:)]);

for xx = 1:Priv.dataSize(nrSession,1)
  if nrBand == 0
    if ~strcmp(en,Priv.computerInfo.endian)
      A = swapbytes(memMaps.origMap.([Priv.prefixSession num2str(nrSession)]...
      ).([Priv.prefixSubject num2str(m)]).Data(xx).tyz);%(:,yy(iter),zz(iter));
    else
      A = memMaps.origMap.([Priv.prefixSession num2str(nrSession)]...
      ).([Priv.prefixSubject num2str(m)]).Data(xx).tyz;%(:,yy(iter),zz(iter));
    end
  else
    if ~strcmp(en,Priv.computerInfo.endian)
      A = swapbytes(memMaps.filtMap.([Priv.prefixSession num2str(nrSession)]...
      ).([Priv.prefixSubjectFilt num2str(m)]).([Priv.prefixFreqBand ...
      num2str(nrBand)]).Data(xx(iter)).tyz;%(:,yy(iter),zz(iter));
    else
      A = memMaps.filtMap.([Priv.prefixSession num2str(nrSession)]...
      ).([Priv.prefixSubjectFilt num2str(m)]).([Priv.prefixFreqBand ...
      num2str(nrBand)]).Data(xx(iter)).tyz;%(:,yy(iter),zz(iter));
    end
  end
  B(xx,:,:,:) = A([2 3 1]);
end
