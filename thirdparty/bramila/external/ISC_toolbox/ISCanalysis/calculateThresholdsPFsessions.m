function calculateThresholdsPFsessions(Params)

% Compute critical thresholds from sumZPF null distributions for across
% session comparison.
%
% See also: RUNANALYSIS, INITPARAMS

% Last updated: 21.6.2011 by Jukka-Pekka Kauppi
% Tampere University of Technology
% Department of Signal Processing
% e-mail: jukka-pekka.kauppi@tut.fi

Priv = Params.PrivateParams;
Pub = Params.PublicParams;
sessComps = ((Priv.nrSessions)^2-(Priv.nrSessions))/2;

critVals = [0.05 0.01 0.001];
r = 0;
for s = 0:Pub.nrFreqBands
    for k = 1:sessComps
        try
            load([Priv.PFsessionDestination 'band' num2str(s) 'valsPFsessComp' num2str(k)])
            vals1 = vals1(:)';
            vals2 = vals2(:)';
            vals = [vals1 vals2];
            vals = fliplr(sort(vals));
            for m = 1:length(critVals)
                Th(m) = vals(1 + floor(length(vals)*critVals(m)));
            end
            save([Priv.PFsessionDestination 'ThPFBand' num2str(s) Priv.prefixSessComp num2str(k) 'win' num2str(r)],'critVals','Th')
        catch
            warning(['Could not find sum ZPF permutation distribution for session difference (in frequency band ' num2str(s) '), comparison ' num2str(k) ' -> significance test skipped.'])
        end
    end
end


% set session contrast maps non-writable:
% load([Pub.dataDestination 'memMaps'])
% for k = 0:Pub.nrFreqBands
%     for m = 1:sessComps
%         try
%             memMaps.(Priv.PFMapSessionName).whole.([Priv.prefixFreqBand num2str(k)]...
%                 ).cor.([Priv.prefixSessComp num2str(m)]).Writable = false;
%             memMaps.(Priv.PFmatMapSessionName).whole.([Priv.prefixFreqBand num2str(k)]...
%                 ).cor.([Priv.prefixSessComp num2str(m)]).Writable = false;
%         catch
%             lasterr
%         end
%     end
% end
% save([Pub.dataDestination 'memMaps'],'memMaps')
% 
