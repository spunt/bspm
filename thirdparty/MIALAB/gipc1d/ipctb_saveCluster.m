function bOk = ipctb_saveCluster(sClustFile, sFilesOrig, sDirMain, nSubj, nSess)
global gIpc;
bOk=false;

%File that contains the voxels of interest. 
%file = [sDirMain gIpc.sFileT2 '.img'];

iCorrelVoxAll=ipctb_spm_read_vols(ipctb_spm_vol(gIpc.sThreshMaskFile)); %it is threshold mask and not structural that is needed here
iCorrelVoxAll(isfinite(iCorrelVoxAll)==0) = 0;
iCorrelVoxAll = find(iCorrelVoxAll ~= 0);

% temp=reshape(temp,gIpc.xdim*gIpc.ydim*gIpc.zdim,1);
% ind = (temp > mean(temp));

% %Create Mask for Voxels of Interest
% rho = ipctb_spm_read_vols(ipctb_spm_vol(file));
% rho(isfinite(rho)==0) = 0;
%iCorrelVoxAll = find(rho>gIpc.dTMin2Sam); %was 3

if ~isempty(iCorrelVoxAll)
    %Load Subject Data and detrend/normalize
    h = waitbar(0,['Processing clusters...']);
    codMeans = zeros(nSubj*nSess,1);
    adCorrelMasked = zeros(size(iCorrelVoxAll, 1), gIpc.tdim*nSess);    
    for i = 1:nSubj;
        for iSess = 1:nSess
            adTemp = zeros(size(iCorrelVoxAll, 1), gIpc.tdim);                        
            for iTime = 1:gIpc.tdim
                % disp(['Loading dataset ',num2str(iTime)]);
                adVolTimeSlice = ipctb_spm_read_vols(ipctb_spm_vol(sFilesOrig(1,iSess+(i-1)*nSess).name(iTime,:)));
                adVolTimeSlice(isfinite(adVolTimeSlice)==0) = 0;
                adVolTimeSlice = reshape(adVolTimeSlice, gIpc.xdim*gIpc.ydim*gIpc.zdim, 1);
                adTemp(:, iTime) = adVolTimeSlice(iCorrelVoxAll);                
            end
            codMeans(iSess + (i-1)*nSess) = mean(mean(adTemp,2));
            adTemp(:,:) = detrend(adTemp')';
            adCorrelMasked(:,(1+(iSess-1)*gIpc.tdim):(iSess*gIpc.tdim)) = adCorrelMasked(:,(1+(iSess-1)*gIpc.tdim):(iSess*gIpc.tdim)) + adTemp;
        end        
        waitbar(i/(nSubj*nSess),h);
    end
    dMeanOrig = mean(codMeans);
    close(h);
    
    adCorrelMasked = 100 * adCorrelMasked / dMeanOrig; %signal change percent relative to the mean

    % Clustering ...
%     [idc, adTcCat] = kmeans(anCorrelMaskedGrp,gIpc.nClusters,'Replicates',5); % The clustering averaged/replicated 5 times
    [idc, adTcCat, dummy] = ipctb_pir_kmeans2(adCorrelMasked,gIpc.nClusters);    
    save([sDirMain sClustFile 'Copy.mat'], 'adCorrelMasked', 'idc');
    clear adCorrelMasked;

    cluster_image = zeros(gIpc.xdim*gIpc.ydim*gIpc.zdim,1);
    cluster_image(iCorrelVoxAll) = idc;
    cluster_vol = reshape(cluster_image,gIpc.xdim,gIpc.ydim,gIpc.zdim);

    %save Time Courses for each cluster, timepoint and session in 3d matrix
    adTc = zeros(gIpc.nClusters, gIpc.tdim, nSess);
    for iSess = 1:nSess
        adTc(:,:,iSess) = adTcCat(:, (1+(iSess-1)*gIpc.tdim):(iSess*gIpc.tdim));
    end
    save([sDirMain sClustFile '_tc.mat'], 'adTc');

    % Save the clustered brain
    gIpc.V.fname=[sDirMain sClustFile '.img'];
    sOldDir = pwd;
    cd(gIpc.sDirSpm);
    ipctb_spm_write_vol(gIpc.V, cluster_vol);
    cd(sOldDir);    
    disp('Successfully saved cluster');
else
    disp('No cluster saved since no voxels were significant!');
end
bOk=true;