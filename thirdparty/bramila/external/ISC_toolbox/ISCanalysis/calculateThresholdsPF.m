function calculateThresholdsPF(Params)

% 
%
%


Priv = Params.PrivateParams;  
Pub = Params.PublicParams;
freqComps = ((Priv.maxScale+2)^2-(Priv.maxScale+2))/2;

% for s = 1:Priv.nrSessions
%     for k = 1:freqComps
%         Fi = [Priv.PFDestination 'session' num2str(s) 'valsPFfreqComp' num2str(k) '.mat'];
%         if exist(Fi) == 2
%             disp(['File ' Fi ' already exists, skipping threshold calculation...'])
%             return
%         end
%     end
% end

critVals = [0.05 0.005 0.001];
r = 0;
for s = 1:Priv.nrSessions
  for k = 1:freqComps
    try
      load([Priv.PFDestination 'session' num2str(s) 'valsPFfreqComp' num2str(k)])
      vals1 = vals1(:)';
      vals2 = vals2(:)';
      vals = [vals1 vals2];
      vals = fliplr(sort(vals));
      for m = 1:length(critVals)
        Th(m) = vals(1 + floor(length(vals)*critVals(m)));
      end
      save([Priv.PFDestination 'ThPFSession' num2str(s) Priv.prefixFreqComp num2str(k) 'win' num2str(r)],'critVals','Th')
    catch err
      warning(err, ['Could not find sum ZPF distribution for session ' num2str(s) ', comparison ' num2str(k) ' -> ignored.'])
    end
  end
end


% set frequency-band contrast maps non-writable:
% load([Pub.dataDestination 'memMaps'])
% for k = 1:Priv.nrSessions
%     for m = 1:freqComps
%         try
%         memMaps.(Priv.PFMapName).whole.([Priv.prefixSession num2str(k)]...
%             ).cor.([Priv.prefixFreqComp num2str(m)]).Writable = false;
%         memMaps.(Priv.PFmatMapName).whole.([Priv.prefixSession num2str(k)]...
%             ).cor.([Priv.prefixFreqComp num2str(m)]).Writable = false;
%         catch err
%             disp(err.message)
%         end
%     end
% end
% save([Pub.dataDestination 'memMaps'],'memMaps')