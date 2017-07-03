function memMaps = byteSwapping(memMaps,Params)

Priv = Params.PrivateParams;
Pub = Params.PublicParams;

disp('Swapping time-series:')
for nrSession = 1:Priv.nrSessions
  disp(['session: ' num2str(nrSession) '/' num2str(Priv.nrSessions)])
  for nrBand = 1:Priv.maxScale + 1
    disp(['band: ' num2str(nrBand) '/' num2str(Priv.maxScale + 1)])
    for k = 1:Priv.nrSubjects
      disp(['subject: ' num2str(k) '/' num2str(Priv.nrSubjects)])
      if nrBand == 0 % load full-band data
        H = memMaps.(Priv.origMapName).([Priv.prefixSession ...
        num2str(nrSession)]).([Priv.prefixSubject ...
        num2str(k)]);
      else % load sub-band data
        H = memMaps.(Priv.filtMapName).([Priv.prefixSession ...
        num2str(nrSession)]).([Priv.prefixSubjectFilt ...
        num2str(k)]).([Priv.prefixFreqBand num2str(nrBand)]);
      end
      H.Writable = logical(1);
      for xx = 1:Priv.dataSize(nrSession,1)
        % if k == 2
        %   if xx >= 30 && xx <= 56 
        %     continue
        %   end
        % end
        size(H)
        whos H
        k
        xx
        nrBand
        cDat = swapbytes(H.Data(xx).tyz(:,:,:));
        H.Data(xx).tyz(:,:,:) = cDat;
        
        if nrBand == 0 % load full-band data
          memMaps.(Priv.origMapName).([Priv.prefixSession ...
          num2str(nrSession)]).([Priv.prefixSubject ...
          num2str(k)]) = cDat;
        else
          memMaps.(Priv.filtMapName).([Priv.prefixSession ...
          num2str(nrSession)]).([Priv.prefixSubjectFilt ...
          num2str(k)]).([Priv.prefixFreqBand num2str(nrBand)]) = H;
        end
        min(cDat(:))
        max(cDat(:))
      end
      H.Writable = logical(0);
      %            disp(['min value after swapping: ' num2str(min(mi))])
      %            disp(['max value after swapping' num2str(max(Ma))])
      %            clear mi Ma
      
    end
  end
end


      % fn = {'whole','win'};
      % 
      % R(1) = Params.PublicParams.ssiOn;
% R(2) = Params.PublicParams.nmiOn;
% R(3) = Params.PublicParams.corOn;
% R(4) = Params.PublicParams.kenOn;
% R = nonzeros(R.*(1:4));
% M = -1e100;
% mi = 1e100;
% 
% disp('Swapping result-maps:')
% 
% for k = 1:length(fn)
%     for m = 0:Params.PrivateParams.maxScale + 1
%         for n = 1:Params.PrivateParams.nrSessions
%             for p = 1:length(R)
%                 H = memMaps.(Params.PrivateParams.resultMapName).(fn{k}).([Params.PrivateParams.prefixFreqBand ...
%                     num2str(m)]).([Params.PrivateParams.prefixSession ...
%                     num2str(n)]).(Params.PrivateParams.simM{R(p)});
%                 H.Writable = logical(1);
%                 disp([num2str(k) ' ' num2str(m) ' ' num2str(n) ' ' num2str(p)])
%                 for d = 1:length(H.Data)
%                     H.Data(d).xyz = swapbytes(H.Data(d).xyz);
%                     Mc = max(max(max(H.Data(d).xyz(:,:,:))));
%                     mc = min(min(min(H.Data(d).xyz(:,:,:))));
%                     if Mc > M
%                         M = Mc;
%                     end
%                     if mc < mi
%                         mi = mc;
%                     end
%                 end
%                 mi
%                 M
%                 H.Writable = logical(0);
%                 memMaps.(Params.PrivateParams.resultMapName).(fn{k}).([Params.PrivateParams.prefixFreqBand ...
%                     num2str(m)]).([Params.PrivateParams.prefixSession ...
%                     num2str(n)]).(Params.PrivateParams.simM{R(p)}) = H;
%             end
%         end
%     end
% end
% 
