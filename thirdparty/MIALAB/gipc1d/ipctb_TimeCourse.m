function bOk = ipctb_TimeCourse(sTitle, sFileRead, sFileWrite)
% draws the time locked time courses (last presentational step)

    global gIpc;
    bOk = false;
    
    % time points per second after intrapolation
    dIntraTpsPerSec = gIpc.dSyncFact/gIpc.dTr;    
        
    % get adTc (3d matrix with clusters, timepoints and sessions)
    load (sFileRead);
    
    % Validation
    if size(adTcLockedMn, 1) > size(gIpc.cmClust,1)
        msgbox(sprintf('Not enough colors in gIpc.cmClust to display time courses (check gipc_defaults.m).\n ---------------------------------------------------------------'));
    end    
    
    nSess = size(adTcLockedMn,3);
    for iSess=1:nSess
        hPlot = figure;
                
        for i = 1:size(adTcLockedMn, 1) % Loop Clusters
            dBaseline = adTcLockedMn(i,1, iSess);
            plot(   (1:size(adTcLockedMn, 2)) / dIntraTpsPerSec, adTcLockedMn(i,:, iSess)-dBaseline, 'LineWidth',1,'Color',gIpc.cmClust(i+1,:)); 
            hold on;
        end
        
    %   legend('1', '2', '3', '4', '5');
        xlabel('Seconds'); 
        ylabel('% BOLD signal change / undetrended mean of BOLD'); 
        title(['Timelock on onset "' gIpc.sOnsLockType '" using ' sTitle ', ' num2str(nSess) ' session(s), ' num2str(gIpc.nClusters) ' cluster(s)']);
        'Sz, N=3, only 1st session, 5 clusters'; 
        axis tight;
        saveas(hPlot, sFileWrite);
        close(hPlot);
    end
    
    bOk = true;