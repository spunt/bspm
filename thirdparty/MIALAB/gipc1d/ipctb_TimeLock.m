function bOk = ipctb_TimeLock(sFileRead, sFileWrite)

    global gIpc;
    bOk = false;
        
    % get adTc (3d matrix with clusters, timepoints and sessions)
    load (sFileRead);
    % Intrapolate timecourses to match double resolution of onsets 
    adTc = ipctb_3dInterp1(adTc, gIpc.dSyncFact, 2);      
    
    % time points per second after intrapolation
    dIntraTpsPerSec = gIpc.dSyncFact/gIpc.dTr;
    
    % Validation
    if size(adTc, 1) > size(gIpc.cmClust,1)
        msgbox(sprintf('Not enough colors in gIpc.cmClust to display time courses (check gipc_defaults.m).\n ---------------------------------------------------------------'));
    end    
    
    nSess = size(adTc,3);

    %number of sessions
    num_sess = size(adTc,3);
    %number of clusters
    num_clusters = gIpc.nClusters;

    % Last possible time point for time lock. Intrapolated timecourse / dIntraTpsPerSec
    dScanEndSec = size(adTc,2) / dIntraTpsPerSec;

    nTotOnsets = 0;
    suTemp(num_sess).rodOnsets = [];
    for iSess = 1:num_sess %loop through each session.      
        for iType = 1:size(gIpc.nSess(iSess).susOnset,1)
            if strcmp(gIpc.nSess(iSess).susOnset{iType,1}, gIpc.sOnsLockType)
                rodTemp = load([gIpc.sDirOnset gIpc.nSess(iSess).susOnset{iType,3}]);
                if gIpc.bOnsetInTr0Sec1 == 0
                    rodTemp = rodTemp/gIpc.dTr; %changes onsets to sec if needed                
                end
                irodTemp = find( rodTemp <= (dScanEndSec - ceil(gIpc.dLockTime)) );
                suTemp(iSess).rodOnsets = [suTemp(iSess).rodOnsets rodTemp(irodTemp)];
            end
        end
        nTotOnsets = nTotOnsets + size(suTemp(iSess).rodOnsets,2);
    end

    nTcLen = floor(gIpc.dLockTime * dIntraTpsPerSec);
    adTcTimeLocked = zeros(nTotOnsets, nTcLen);
    adTcLockedMn = zeros(num_clusters, nTcLen);
    iTc = 1;
    for iClust = 1:num_clusters,
        for iSess = 1:num_sess, %loop through each session.
            for iLocks = 1:size(suTemp(iSess).rodOnsets,2),
                % convert time to index (rounded to closest index)
                iTimeLock = round((suTemp(iSess).rodOnsets(iLocks)) * dIntraTpsPerSec);
                adTcTimeLocked(iTc, :) = adTc(iClust, iTimeLock:(iTimeLock+nTcLen - 1), iSess);
                iTc = iTc + 1;
            end
        end
        adTcLockedMn(iClust,:) = mean(adTcTimeLocked, 1);
    end
    
    save(sFileWrite, 'adTcLockedMn');
    bOk = true;