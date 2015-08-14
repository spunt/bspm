function bOk = ipctb_saveClusterTcConverse(sClustFile1, sClustFile2, sDirMain, nSess)
    global gIpc;
    bOk = false;

    load([sDirMain sClustFile1 'Copy.mat'], 'idc'); %cluster mask
    load([sDirMain sClustFile2 'Copy.mat'], 'adCorrelMasked'); % time courses for each voxel
    % 'adCorrelMasked', 'idc', 'iCorrelVoxAll'

    adTcCat = zeros(gIpc.nClusters, gIpc.tdim * nSess);
    for iClust = 1:gIpc.nClusters
        iMask2 = find(idc == iClust);
        adTcCat(iClust, :) = mean(adCorrelMasked(iMask2,:));
    end

    %save Time Courses for each cluster, timepoint and session in 3d matrix
    adTc = zeros(gIpc.nClusters, gIpc.tdim, nSess);
    for iSess = 1:nSess
        adTc(:,:,iSess) = adTcCat(:, (1+(iSess-1)*gIpc.tdim):(iSess*gIpc.tdim));
    end
    save([sDirMain sClustFile1 '_tcConverse.mat'], 'adTc');

    bOk = true;