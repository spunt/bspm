function bOk = ipctb_saveCorrMaps(sPathExport, suFiles, nSubjTot, nSess)
% function that loads specified img-files and 
global gIpc;
bOk = false;

k=0;
% nTotVol = gIpc.xdim*gIpc.ydim*gIpc.zdim;
nTotVol = gIpc.xdim*gIpc.ydim;

% nSubjTot are extracted from srinivas file getter for group1 and
% group2. Here is also the session information. That should be that
% Can I move out this function to separate file?

savefilenames=[]; %Create an empty array for filling.
%Generate the list of filenames for saving your new data. They will be in
%the form rho(subject#)vs(subject#) ie. rho01vs02. 
for i = 1:(nSubjTot);
    if (i<10);
        new_i=['0',num2str(i)];
    else
        new_i=[num2str(i)];
    end;   
    for j = (i+1):(nSubjTot);
        if (j<10);     
            new_j=['0',num2str(j)];
        else
            new_j=[num2str(j)];
        end;
    tmp = [ipctb_backslash(sPathExport),'rho',num2str(new_i),'vs',num2str(new_j)];
    savefilenames = cat(1,tmp,savefilenames);
    end;
end;

%We're creating a threshold here using the first subject's data and
%creating a mask that will be used to filter for in-brain voxels only.
disp(['Creating thresholding mask']);
%temp=spm_read_vols(spm_vol('/export/research/analysis/human/collaboration/fbirn/users/dkim/ipc_aod/sampleEPI.img'));
temp=ipctb_spm_read_vols(ipctb_spm_vol(gIpc.sThreshMaskFile));
temp(isfinite(temp)==0) = 0;
temp=reshape(temp,gIpc.xdim*gIpc.ydim*gIpc.zdim,1);
ind = (temp > mean(temp));

subj1=zeros(nTotVol * nSess, gIpc.tdim);
subj2=zeros(nTotVol * nSess, gIpc.tdim);

for i = 1:nSubjTot;
    %First subject and first run.
        
    for j = (i+1):nSubjTot;
        %Second subject and first run.
        rho = zeros(gIpc.xdim*gIpc.ydim*gIpc.zdim, 1);
        
        
        h = waitbar(0,['Processing subject ' num2str(i) ' (of ' num2str(nSubjTot-1) ') combinations...']);
        for iSlice = 1:gIpc.zdim

            for iSess = 1:nSess
%                 disp(['Loading subject ',num2str(i),' run ',num2str(iSess)]);       
                VIpc = ipctb_spm_vol(suFiles(1,iSess+(i-1)*nSess).name);
                for iLp = 1:gIpc.tdim
                    try
                        temp2 = reshape(ipctb_spm_slice_vol(VIpc(iLp), ipctb_spm_matrix([0 0 iSlice]), [gIpc.xdim gIpc.ydim], 0) , nTotVol , 1);
                        temp2(isfinite(temp2)==0) = 0;
                        subj1(1+(iSess-1)*nTotVol:iSess*nTotVol, iLp)=temp2;
                    catch
                        msgbox(sprintf(['Your image time series does not match\nbetween the two groups.\n ---------------------------------------------------------------']));
                        return;                        
                    end
                end
                clear temp2;
                subj1(1+(iSess-1)*nTotVol:iSess*nTotVol, :)=detrend(subj1(1+(iSess-1)*nTotVol:iSess*nTotVol, :)')';
                subj1(1+(iSess-1)*nTotVol:iSess*nTotVol, :)=subj1(1+(iSess-1)*nTotVol:iSess*nTotVol, :)./repmat(std(subj1(1+(iSess-1)*nTotVol:iSess*nTotVol, :),[],2),1, gIpc.tdim);
                
            end    

            %Second subject     
            for iSess = 1:nSess
                % disp(['Loading subject ',num2str(j),' run ',num2str(iSess)]);
                VIpc = ipctb_spm_vol(suFiles(1,iSess+(j-1)*nSess).name);
                for iLp = 1:gIpc.tdim
                    temp2 = reshape(ipctb_spm_slice_vol(VIpc(iLp), ipctb_spm_matrix([0 0 iSlice]), [gIpc.xdim gIpc.ydim], 0) , nTotVol , 1);
                    temp2(isfinite(temp2)==0) = 0;
                    subj2(1+(iSess-1)*nTotVol:iSess*nTotVol, iLp)= temp2;
                    clear temp2;
                end
                subj2(1+(iSess-1)*nTotVol:iSess*nTotVol, :)=detrend(subj2(1+(iSess-1)*nTotVol:iSess*nTotVol, :)')';
                subj2(1+(iSess-1)*nTotVol:iSess*nTotVol, :)=subj2(1+(iSess-1)*nTotVol:iSess*nTotVol, :)./repmat(std(subj2(1+(iSess-1)*nTotVol:iSess*nTotVol, :),[],2),1, gIpc.tdim);
                
            end    

            % disp(['Regressing subject ',num2str(i),' with subject ' num2str(j)]);
            %The for loop runs through the entire voxel space of a single
            %run for a single subject within the mask created by in. 
            %In other words, work only within the brain and not waste time on
            %other voxels.        
            aIdX = [zeros(nSess*gIpc.tdim, nSess) ones(nSess*gIpc.tdim,1)];
            adSubj2 = zeros(nSess*gIpc.tdim,1);
            for x = 1 : nTotVol
                if (rem(x, 10000) == 0), disp([num2str(x) ' out of ' num2str(nTotVol)]);
                end;
                if ((ind(    (x   + (iSlice-1)*nTotVol ),1) ~= 0)   ); % && min(x ~= 3035:3041)
                    %We are going to build a design matrix to allow for variances
                    %in the run which we might not have accounted for if we just
                    %concatenated the two runs.
                    %Create the Design Matrix.
                    for iSess = 1:nSess
                        aIdX( 1+(iSess-1)*gIpc.tdim:iSess*gIpc.tdim , iSess) = subj1(x+(iSess-1)*nTotVol, :)';
                        adSubj2( 1+(iSess-1)*gIpc.tdim:iSess*gIpc.tdim,1) = subj2(x+(iSess-1)*nTotVol,:)';
                    end
                    %Perform the Regression.
                    [b, stats]=ipctb_ica_regress(adSubj2,aIdX);
                    rho(x+(iSlice-1)*nTotVol,1) = ((sign(b(1)+b(2))) * sqrt(abs(stats(1))));
                end;
            end

            waitbar(iSlice/gIpc.zdim,h)
        
        end
        close(h)

        rho = reshape(rho, gIpc.xdim, gIpc.ydim, gIpc.zdim); %Reshape the correlation map.
        %This part displays the status of the correlation script.
        k = k+1;
        disp('  ');
        disp([num2str(k) ' FILE(S) WRITTEN OUT OF ' num2str(size(savefilenames),1)]);
        %savefilenames matrix
        save (savefilenames(k, :),  'rho');
        %clears some variables so as to not to get the dreaded "out of
        %memory" message
        clear subj2;
        clear subj2;
    end;
    clear subj1;
    clear subj1;
end;
bOk = true;