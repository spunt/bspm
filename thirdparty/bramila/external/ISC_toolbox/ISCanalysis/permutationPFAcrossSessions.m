function permutationPFAcrossSessions(Params,nrBand,sessBlock)
% This function generates null distribution for sum ZPF
% statistic. The sign of the ZPF statistic of the subject
% pairs is randomly flipped and summed up several times. At
% single iteration, same randomization is used across all brain
% voxels and maximum across the voxels is selected to obtain final
% realization. This procedure takes automatically care of the multiple
% comparisons.
%
% See also: RUNANALYSIS, INITPARAMS

% Last updated: 21.6.2011 by Jukka-Pekka Kauppi
% Tampere University of Technology
% Department of Signal Processing
% e-mail: jukka-pekka.kauppi@tut.fi


showTime(1);

Priv = Params.PrivateParams;
Pub = Params.PublicParams;



if Pub.sessionCompOn == 0
    disp('Sum ZPF map computation across sessions not selected...')
    return
end

if ~Pub.corOn
    disp('Correlation coefficient needs to be selected in order to compute sum ZPF maps...')
    return
end

if ~Pub.calcStandard
    disp('Correlation coefficient needs to be selected in order to compute sum ZPF maps...')
    return
end

for xx = 1:Priv.nrSessions
    za(xx) = Priv.dataSize(xx,4);
end
if length(unique(za)) ~= 1
    disp('Each session must have same number of data points in order to compute sum ZPF maps across sessions, calculation skipped...')
    return
end

if Priv.nrSessions < 2
    disp('Only one session specified, sumZPF across sessions cannot be computed...')
    return
end

for fc = 1:length(sessBlock)
    Fi = [Priv.PFsessionDestination 'band' num2str(nrBand) 'valsPFsessComp' num2str(sessBlock(fc)) '.mat'];
    if exist(Fi) == 2
        disp(['File ' Fi ' already exists, sum ZPF permutation test skipped...'])
        return
    end
end



load([Pub.dataDestination 'memMaps'])

nrPerm = Pub.permutSessionComp;
nrSubjectPairs = (Priv.nrSubjects^2-Priv.nrSubjects)/2;

% load mask:
if Pub.useTemplate
  bmask = load_nii(Priv.brainMask);
  bmask = single(bmask.img);
else
  bmask = load(Priv.brainMask);
  fiel = fields(bmask);
  bmask = bmask.(fiel{1});
  bmask = single(bmask);
end
bmask = logical(bmask);

% remove zero-planes from data to easy computational burden:

% find zero-planes:
IDX = [1 1 1 1 1 1];
for s = 1:Priv.dataSize(1,1)
    if sum(sum(squeeze(bmask(s,:,:)))) == 0
        IDX(1) = IDX(1)+1;
    else
        break
    end
end
for s = Priv.dataSize(1,1):-1:1
    if sum(sum(squeeze(bmask(s,:,:)))) == 0
        IDX(2) = IDX(2)+1;
    else
        break
    end
end
for s = 1:Priv.dataSize(1,2)
    if sum(sum(squeeze(bmask(:,s,:)))) == 0
        IDX(3) = IDX(3)+1;
    else
        break
    end
end
for s = Priv.dataSize(1,2):-1:1
    if sum(sum(squeeze(bmask(:,s,:)))) == 0
        IDX(4) = IDX(4)+1;
    else
        break
    end
end
for s = 1:Priv.dataSize(1,3)
    if sum(sum(squeeze(bmask(:,:,s)))) == 0
        IDX(5) = IDX(5)+1;
    else
        break
    end
end
for s = Priv.dataSize(1,3):-1:1
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

% do permutation for each session comparison pair given in sessionComp:
for fc = 1:length(sessBlock)

    disp(['Session comparison ' num2str(sessBlock(fc)) '....'])

    % get ZPF data matrices for each subject pair:
    PFdata = memMaps.(Priv.PFmatMapSessionName).whole.([...
        Priv.prefixFreqBand num2str(nrBand)]).cor.([...
        Priv.prefixSessComp num2str(sessBlock(fc))]).Data.xyzc;    

    PFsmaller = PFdata(IDX(1):end-(IDX(2)-1),IDX(3):end-(IDX(4)-1),...
        IDX(5):end-(IDX(6)-1),:);

    siz = prod(size(bmask));
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
    % save null distributions of minimal and maximal sumZPF statistic:
    save([Priv.PFsessionDestination 'band' num2str(nrBand) 'valsPFsessComp' num2str(sessBlock(fc))],'vals1','vals2')
end

showTime(0);
