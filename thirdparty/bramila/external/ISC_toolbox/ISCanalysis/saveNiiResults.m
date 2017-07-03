function saveNiiResults(Params)
% function generates nii files from the computed results
%

Priv = Params.PrivateParams;
Pub = Params.PublicParams;

%load Memmaps and one original header to get the voxel dimensions
load([Pub.dataDestination 'memMaps'])
%nii_orighdr = load_nii_hdr(Pub.subjectSource{1});

for nrBand = 0:Pub.nrFreqBands
  for nrSession = 1:Priv.nrSessions
        %read the data from the memmaps
        img=memMaps.resultMap.whole.(['band' num2str(nrBand)]).(['Session' num2str(nrSession)]).cor.Data.xyz;

        % generate nii structure
        %nii_out = make_nii(img, nii_orighdr.dime.pixdim(2:4),[],[],'ISCtoolbox');
        nii_out = make_nii(img, Priv.dataSize(nrSession,1:3) ,[],[],'ISCtoolbox'); 

        % save the nii to the result folder of destination
        save_nii(nii_out,[Priv.resultsDestination 'ISCcorrmapBand' num2str(nrBand) 'Session' num2str(nrSession) '.nii']);
  end
end
clear memMaps;