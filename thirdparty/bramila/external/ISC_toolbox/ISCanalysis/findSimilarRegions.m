function [r05,r005,R05,R005] = findSimilarRegions(Params)

load([Params.PublicParams.dataDestination 'memMaps'])
Priv = Params.PrivateParams;
Pub = Params.PublicParams;
R005 = cell(Priv.nrSessions,Priv.maxScale+2,max(Priv.nrTimeIntervals));
R05 = cell(Priv.nrSessions,Priv.maxScale+2,max(Priv.nrTimeIntervals));
r005 = cell(Priv.nrSessions,Priv.maxScale+2);
r05 = cell(Priv.nrSessions,Priv.maxScale+2);
for k = 0:Priv.maxScale+1
  for m = 1:Priv.nrSessions
    disp(['Session ' num2str(m) ', freq.band ' num2str(k)])
    load([Priv.resultsDestination 'Th' Priv.prefixFreqBand... 
    num2str(k) Priv.prefixSession num2str(m) 'win0'],'Th')
    q = memMaps.resultMap.whole.([Priv.prefixFreqBand num2str(k)]...
          ).([Priv.prefixSession num2str(m)]).cor.Data.xyz;
    r05{m,k+1} = find(q >= Th(2));
    r005{m,k+1} = find(q >= Th(6));
    load([Priv.resultsDestination 'Th' Priv.prefixFreqBand... 
    num2str(k) Priv.prefixSession num2str(m) 'win1'],'Th')
    for t = 1:Priv.nrTimeIntervals(m)
      Q = memMaps.resultMap.win.([Priv.prefixFreqBand num2str(k)]...
      ).([Priv.prefixSession num2str(m)]).cor.Data(t).xyz;
      R05{m,k+1,t} = find(Q >= Th(2));
      R005{m,k+1,t} = find(Q >= Th(6));
    end    
  end
end
save([Priv.resultsDestination 'Regions'],'R005','R05','r005','r05')
