function permutationPF(Params,nrSession,freqComp)
% This function generates null distribution for sum ZPF
% statistic. The sign of the ZPF statistic of the subject
% pairs is randomly flipped and summed up several times. At
% single iteration, same randomization is used across all brain
% voxels and maximum across the voxels is selected for final
% realization. This procedure takes automatically care of the multiple % comparisons.
%
% See also: RUNANALYSIS, INITPARAMS

% Last updated: 27.1.2010 by Jukka-Pekka Kauppi
% Tampere University of Technology
% Department of Signal Processing
% e-mail: jukka-pekka.kauppi@tut.fi


showTime(1);

Priv = Params.PrivateParams;
Pub = Params.PublicParams;

if ~Pub.corOn
    disp('Correlation criterion not selected, sum Z-Pearson-Filon permutation test skipped...')
    return
end

if ~Pub.calcStandard
    disp('Standard inter-subject synchronization analysis not specified, sum Z-Pearson-Filon permutation test skipped...')
    return
end

if Pub.nrFreqBands == 0
    disp('No frequency subbands selected, sum Z-Pearson-Filon permutation test skipped...')
    return
end

for fc = 1:length(freqComp)
    Fi = [Priv.PFDestination 'session' num2str(nrSession) 'valsPFfreqComp' num2str(freqComp(fc)) '.mat'];
    if exist(Fi,'file') == 2
        disp(['File ' Fi ' already exists, sum ZPF permutation test skipped...'])
        return
    end
end

load([Pub.dataDestination 'memMaps'])

nrPerm = Pub.permutFreqComp;
nrSubjectPairs = (Priv.nrSubjects^2-Priv.nrSubjects)/2;
if strcmp(Pub.fileFormat,'nii')
    bmask = load_nii(Priv.brainMask);
    bmask = single(bmask.img);
elseif strcmp(Pub.fileFormat,'mat')
    bmask = load(Priv.brainMask);
    fiel = fields(bmask);
    bmask = bmask.(fiel{1});
    bmask = single(bmask);
else
    error('Mask must be mat- or nii-file!')
end
bmask = logical(bmask);

% remove zero-planes from data to easy computational burden:

% find zero-planes:
IDX = [1 1 1 1 1 1];
for s = 1:Priv.dataSize(nrSession,1)
    if sum(sum(squeeze(bmask(s,:,:)))) == 0
        IDX(1) = IDX(1)+1;
    else
        break
    end
end
for s = Priv.dataSize(nrSession,1):-1:1
    if sum(sum(squeeze(bmask(s,:,:)))) == 0
        IDX(2) = IDX(2)+1;
    else
        break
    end
end
for s = 1:Priv.dataSize(nrSession,2)
    if sum(sum(squeeze(bmask(:,s,:)))) == 0
        IDX(3) = IDX(3)+1;
    else
        break
    end
end
for s = Priv.dataSize(nrSession,2):-1:1
    if sum(sum(squeeze(bmask(:,s,:)))) == 0
        IDX(4) = IDX(4)+1;
    else
        break
    end
end
for s = 1:Priv.dataSize(nrSession,3)
    if sum(sum(squeeze(bmask(:,:,s)))) == 0
        IDX(5) = IDX(5)+1;
    else
        break
    end
end
for s = Priv.dataSize(nrSession,3):-1:1
    if sum(sum(squeeze(bmask(:,:,s)))) == 0
        IDX(6) = IDX(6)+1;
    else
        break
    end
end
% remove planes:
bmask = bmask(IDX(1):end-(IDX(2)-1),IDX(3):end-(IDX(4)-1),...
    IDX(5):end-(IDX(6)-1));
% vectorize mask:
bmask = bmask(:);

% do permutation for each frequency comparison pair given in freqComp:
for fc = 1:length(freqComp)
    
    disp(['Frequency comparison ' num2str(freqComp(fc)) '....'])
    
    % get ZPF data matrices for each subject pair:
    PFdata = memMaps.(Priv.PFmatMapName).whole.([Priv.prefixSession ...
        num2str(nrSession)]).cor.([Priv.prefixFreqComp num2str(freqComp(fc))]...
        ).Data.xyzc;
    
    PFsmaller = PFdata(IDX(1):end-(IDX(2)-1),IDX(3):end-(IDX(4)-1),...
        IDX(5):end-(IDX(6)-1),:);
    
    siz = numel(bmask);
    PFsmaller = reshape(PFsmaller,siz,nrSubjectPairs);
    PFsmaller(bmask==0,:) = [];
    siz = sum(bmask);
    
    disp('Generating randomization matrix...')
    
    % flip ZPF values randomly (-1 or 1):
    randVals = 2*(rand(nrPerm-1,nrSubjectPairs) >= 0.5)-1;
    % add also the original combination:
    randVals = [ones(1,nrSubjectPairs);randVals];
    
    clear PFdata
    
    vals1 = single(NaN*zeros(1,nrPerm));
    vals2 = vals1;
    reportIntVal = max(1,(round((nrPerm/100))));
    
    disp('Generating null distribution for z-PF matrices...')
    tic
    for k = 1:nrPerm
        r = repmat(randVals(k,:),siz,1);
        Vs = sum(PFsmaller.*r,2);
        vals1(k) = single(max(Vs));
        vals2(k) = single(min(Vs));
        if mod(k,reportIntVal) == 0
            disp(['Iteration: ' num2str(k) '/' num2str(nrPerm)])
            toc
            tic
        end
    end
    toc
    save([Priv.PFDestination 'session' num2str(nrSession) 'valsPFfreqComp' num2str(freqComp(fc))],'vals1','vals2')
end

showTime(0);
