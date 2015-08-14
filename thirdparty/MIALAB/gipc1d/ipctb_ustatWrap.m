function [codRbAll codSe2All, dTMax] = ipctb_ustatWrap(sPathExport, sGrp, sGrpFile, nSubjTot)
global gIpc;

sOldDir = pwd;

nVox = gIpc.xdim*gIpc.ydim*gIpc.zdim;

codGrp = zeros(nVox,1);
codRbAll = zeros(nVox,1);
codSe2All = zeros(nVox,1);
adPartMapsAll = zeros(nVox/gIpc.nSplit,(nSubjTot^2-nSubjTot)/2);

%create matrix containing braincorrelations in columns
disp(' ');
disp('Loading AOD_3T task pair-wise correlation .mat files...');

for iSplit = 1:gIpc.nSplit
    nMap = 1;
    for p = 1 : nSubjTot;
        for q = (p + 1) : nSubjTot;
            disp(['Loading rho' sprintf('%0.2d' ,p) 'vs' sprintf('%0.2d', q) '.mat']);
            eval(['load ' [ipctb_backslash(sPathExport) ipctb_backslash(sGrp)] 'rho' sprintf('%0.2d' ,p) 'vs' sprintf('%0.2d', q) '.mat']);
%             eval(['load /export/research/analysis/human/collaboration/olin/mialab/users/dkim/projects/old_projects/ipc_aod/healthy_ipc/rho' sprintf('%0.2d' ,p) 'vs' sprintf('%0.2d', q) '.mat']);
            coTemp = reshape(rho, gIpc.xdim * gIpc.ydim * gIpc.zdim, 1);
            adPartMapsAll(:,nMap) = coTemp(1+(iSplit-1)*nVox/gIpc.nSplit:iSplit*nVox/gIpc.nSplit,1);
            nMap = nMap + 1;
        end;
    end;
    [codRb codSe2] = ipctb_ustat(adPartMapsAll, [sGrp ' split ' num2str(iSplit) ' (of ' num2str(gIpc.nSplit) ')']);
    codRbAll(1+(iSplit-1)*nVox/gIpc.nSplit:iSplit*nVox/gIpc.nSplit,1) = codRb;
    codSe2All(1+(iSplit-1)*nVox/gIpc.nSplit:iSplit*nVox/gIpc.nSplit,1) = codSe2; 
    codGrp(1+(iSplit-1)*nVox/gIpc.nSplit:iSplit*nVox/gIpc.nSplit,1) = codRb./sqrt(codSe2);
end

%Mask for all brain data needed to discriminate nans
temp=ipctb_spm_read_vols(ipctb_spm_vol(gIpc.sThreshMaskFile)); 
temp(isfinite(temp)==0) = 0;
temp=reshape(temp,gIpc.xdim*gIpc.ydim*gIpc.zdim,1);
ind = find(temp > mean(temp));

% fdr correction
[nSignif, iFdrMask] = ipctb_fdr(2*(1-ipctb_spm_Tcdf(codGrp(ind), nSubjTot-1)), gIpc.dFdrQ1Samp);
codFdrCorrected = zeros(nVox, 1);
if gIpc.bFdr
    if nSignif ~= 0
        codFdrCorrected(ind(iFdrMask), 1) = codGrp(ind(iFdrMask), 1);
    else
        disp('warning - no significant vocels in ipctb_ustatWrap');
    end
end

dTMax = max(codGrp);

% one dimensional t-values to three dimensions in accordance to brain
codGrp=reshape(codGrp,gIpc.xdim, gIpc.ydim, gIpc.zdim);

%save matrix
save([ipctb_backslash(sPathExport) sGrpFile 'Extra.mat'], 'codGrp');
disp(['Successfully saved ' sGrpFile]);

gIpc.V.fname=[ipctb_backslash(sPathExport) sGrpFile '.img'];
cd(gIpc.sDirSpm);
ipctb_spm_write_vol(gIpc.V, codGrp);
cd(sOldDir);
