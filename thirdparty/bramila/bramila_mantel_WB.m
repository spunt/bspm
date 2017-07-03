function results=bramila_mantel_WB(cfg)
% results=bramila_mantel_WB(cfg)
%
% Whole brain Mantel test for similarity matrices. Input 4D nifti file
% needs to contain the upper triangle elements of the similarity matrix
% (same format as the output of ISC toolbox with option "store
% correaltion matrices"). If number of permutations is set to 0, it
% returns only the mantel correlation volume, otherwise it estimates
% p-values for N sampled voxels. Multiple comparisons correction can be
% done using BHFDR or permutation based cluster correction (add
% citation).
% Input parameters: "type" can only be 'pearson' or 'spearman'
%
%   cgf.infile= path to 4D volume, every volume is a pair of subjects' similarity
%   cfg.mask = path to 3D brain mask
% 	cfg.model = similarity matrix used for model
%   cfg.modelNI = models of no interest (to control for confounds)
% 	cfg.type = 'pearson' or 'spearman' (default 'spearman', slower but less prone to outliers)
% 	cfg.iter = number of permutations (default 5000)
%   cfg.cdfN = samples to estimate CDF (default 100)
%	cfg.cdfMethod= 'ksdensity' 'empirical' 'semiparametric' (default empirical)
%   cfg.BHFDR = if you want to use BH FDR, set a q value e.g. 0.05 (if set to 0, it does not perform BHFDR which is the default option)
%   cfg.CFT = cluster forming threshold for cluster correction in p values (default
%   0.05, if set to zero it does not perform cluster correction)
%
% Output:
%
%   results.r= whole brain mantel r values (no correction)
%   results.p= whole brain mantel p values (no correction)
%   results.rcorr = corrected r values
%   results.cdfp;
%   results.cdfr;
% (c) Enrico Glerean 2016 - Brain and Mind Laboratory Aalto University http://becs.aalto.fi/bml/


%% checking input
% add check that we have readable infile

% add check that we have valid mask, matched with infile size

% add check that model fits dimensions of infile

% initialize modelNI to default
if(~isfield(cfg,'modelNI'))
	cfg.modelNI=0;
else
	% add a check that size is matched
end

% initialize type
if(~isfield(cfg,'type'))
	cfg.type='spearman';
end
% check type
if(strcmp(cfg.type,'pearson')==0 && strcmp(cfg.type,'spearman')==0)
    error('types allowed are only pearson and spearman');
end

% initialize iter
if(~isfield(cfg,'iter'))
	cfg.iter=5000;
end

% initialize cdfN
if(~isfield(cfg,'cdfN'))
	cfg.cdfN=100;
end

% initialize cdfMethod
if(~isfield(cfg,'cdfMethod'))
	cfg.cdfMethod='empirical';
end

% initialize BHFDR
if(~isfield(cfg,'BHFDR'))
	cfg.BHFDR=0;
end

% initialise CFT
if(~isfield(cfg,'CFT'))
	cfg.CFT=0.05;
end

%% starting
modelMat=cfg.model;
kk1=size(modelMat);
Nsubj = kk1(1);

ids=find(triu(ones(kk1(1)),1));

if(size(cfg.modelNI,1)==size(modelMat))
    % cleaning the model
    % add if there's more than one models to regress
    [betas ints residu]=regress(modelMat(ids),[cfg.modelNI(ids) ones(size(ids))]);
    tempmodelMat=zeros(size(modelMat));
    tempmodelMat(ids)=residu/max(residu);
    tempmodelMat=tempmodelMat+tempmodelMat'+eye(size(tempmodelMat));
    modelMat=tempmodelMat;
end


% load data
disp(['Loading ' cfg.infile]);
nii=load_nii(cfg.infile);
mask=load_nii(cfg.mask);

% experimental - use a set of voxels with low ISC to compute CDF
tempI=mean(nii.img,4);
maskI=double(tempI<0.0001);
maskI=double(mask.img).*maskI;

% test that they are square
if(kk1(1) ~= kk1(2))	error('model matrix is not square'); end
if( length(ids) ~= size(nii.img,4))  error('model matrix and brain data do not have the same size'); end


% test that they are symmetrical matrices
temp=modelMat -modelMat';
if(sum(temp(:))~=0) error('Model matrix is not symmetrical'); end

% test that they both are similarity or dissimilarity matrices
diag1=sum(diag(modelMat));
if(diag1~=kk1(1)) error('Model matrix is not a similarity matrix'); end



inmask=find(mask.img>0);
data=reshape(nii.img,[],length(ids));
data=data';

if(size(cfg.modelNI,1)==size(modelMat))
    disp('Cleaning the ISC data');
    % cleaning the data
    for vox=1:size(data,2)
        if(~ismember(vox,inmask)) continue; end
        [betas ints residu]=regress(data(:,vox),[cfg.modelNI(find(triu(true(Nsubj),1))) ones(size(find(triu(true(Nsubj),1))))]);
        tempiscmat=zeros(Nsubj);
        data(:,vox)=residu;
        
        
    end
end

%save -v7.3 debug
rout=zeros(size(nii.img,1),size(nii.img,2),size(nii.img,3));
rout(inmask)=corr(data(:,inmask),modelMat(ids),'type',cfg.type);
mask.nii(find(isnan(rout)))=0; % useful for later
inmask=find(mask.img>0);
rout(find(isnan(rout)))=0;

results.r = rout;

save_nii(make_nii(results.r),'debug.nii')
if(~isfield(cfg,'rCFT'))
%% estimate the uncorrected p values
uu=unique(rout(find(maskI>0)));

if(0)
    step=round(length(uu)/(cfg.cdfN/2))
    usampleA=sort(uu(end:-step:1));
    usampleB=sort(uu(end:-1:end-cfg.cdfN/2));
    usampleB=sort(uu(end:-50:1));
    usample=sort(unique([usampleA; usampleB(end:-1:end-cfg.cdfN/2); uu(end:-1:end-10)]));
    
    
    size(usample)
    
else
    %uu(find(uu<=0))=[];
    step=round(length(uu)/(cfg.cdfN));
    usampleA=sort(uu(end:-step:1));
    usample=sort(unique([usampleA]));
end

steps=(1:cfg.cdfN)/cfg.cdfN;
steps=steps.^0.5;
steps=round(steps*length(uu));
uuL=(1:length(uu))/length(uu);
mU=min(uu);
MU=max(uu);
uuL=uuL*(MU-mU)+mU;


usampleL=uuL([1 steps]);
% find best match
for uID=1:length(usampleL)
    tempdist=abs(uu-usampleL(uID));
    win=find(min(tempdist)==tempdist);
    usample(uID,1)=uu(win(1));
end

save debug steps usample

% now we have replaced inmask with maskI, we need to fix also that rout
% should be masked so that only maskI voxels are used to estimate the CDF
% for mantel

allsurro=zeros(cfg.iter,length(usample));
debugvol=zeros(91,109,91,2);
for uID=1:length(usample);
    val=usample(uID);
    disp([num2str(uID) ' ' num2str(val)])
    ids=find(rout==val);
    ids=ids(1);
    [tempx tempy tempz]=ind2sub(size(rout),ids);
    tempiscs=double(squeeze(nii.img(tempx,tempy,tempz,:)));
    
    
    if(size(cfg.modelNI,1)==size(modelMat))
        % cleaning the model
        % add if there's more than one models to regress
        [betas ints residu]=regress(tempiscs,[cfg.modelNI(find(triu(true(Nsubj),1))) ones(size(find(triu(true(Nsubj),1))))]);
        tempiscmat=zeros(Nsubj);
        tempiscmat(find(triu(true(Nsubj),1)))=residu/max(residu);
        tempiscmat=tempiscmat+tempiscmat'+eye(size(tempiscmat));
        
    else
        tempiscmat=zeros(Nsubj);
        tempiscmat(find(triu(true(Nsubj),1)))=tempiscs;
        tempiscmat=tempiscmat+tempiscmat' +eye(Nsubj);
    end
    [r_mantel p_mantel surro]=bramila_mantel_SP(tempiscmat,modelMat,cfg.iter,cfg.type);
    allsurro(:,uID)=surro;
    mmm=mean(tempiscs);
    disp([num2str(p_mantel) ' ' num2str(mmm) ])
    cdf_sample(uID,1)=1-p_mantel;
    results.cdfr(uID,1)= r_mantel;
    %save(['temp/debug' num2str(uID) '.mat'],'tempiscmat','modelMat','tempx','tempy','tempz','mmm') 
    debugvol(tempx,tempy,tempz,1)=r_mantel;
    debugvol(tempx,tempy,tempz,2)=-log10(p_mantel);
end
results.debugvol=debugvol;
% interpolate CDF over the brain
temp=rout(inmask);

results.cdfp=cdf_sample;
% using KSR
if(0)
    %cdf_sample=medfilt1(cdf_sample,11);
    %[cdf_f,cdf_xi] = ksdensity(cdf_sample,usample,'support','positive','function','cdf');%,'bandwidth',0.05);
    
    
    ksre=ksr(results.cdfr,results.cdfp,0.005,cfg.cdfN); % BW as 0.005
    ksre
    
    %all_cdf=interp1([-1; usample; 1],[0; cdf_sample; 1],temp,'linear');
    all_cdf=interp1([-1; ksre.x'; 1],[0; ksre.f'; 1],temp,'linear');
    results.ksre=ksre;
end

% using empirical CDF
if(0)
    stepfun=zeros(size(results.cdfp));

    meanpoint=length(find(results.cdfp<0.5));
    for tt=1:length(stepfun)
        if(tt<=meanpoint)
            stepfun(tt)=min(results.cdfp(1:tt)); 
        else
            stepfun(tt)=max(results.cdfp(1:tt)); 
        end
    end
    results.stepfun=stepfun;


    tempR=results.cdfr;
    tempR(find(diff(results.cdfr)==0))=[];
    tempstepfun=stepfun;
    tempstepfun(find(diff(results.cdfr)==0))=[];
    all_cdf=interp1([-1; tempR; 1],[0; tempstepfun; 1],temp,'linear');
end

% using semi parametric pareto from permuted mantels
if(1)
    pfit = paretotails(allsurro(:),0.2,0.8);
    all_cdf=cdf(pfit,temp); 
    results.pfit=pfit;
end


end

%% BHFDR
if(cfg.BHFDR>0)
    q=mafdr((1-all_cdf),'BHFDR','true');
    qvol=zeros([size(nii.img,1),size(nii.img,2),size(nii.img,3)]);
    pvol=qvol;
    qvol(inmask)=1-q;
    pvol(inmask)=all_cdf;
    results.q=qvol;
    results.p=pvol;
    % add here code for volume at q< cfg.BHFDR
end

%% cluster correction
ids=find(triu(ones(kk1(1)),1));

if(cfg.CFT>0)
    if(isfield(cfg,'rCFT'))
        rCFT=cfg.rCFT;
    else
        rCFT=min(temp((find((1-all_cdf)<cfg.CFT)))) % cluster forming threshold at cfg.CFT
    end
    % generate surrogate models, first iter is untouched
    surro=zeros(length(ids),cfg.CCiter);
    for i = 1 :cfg.CCiter
        shuf=randperm(Nsubj);
        if(i==1) shuf=1:Nsubj; end
        temp=modelMat(shuf,shuf);
        surro(:,i)=temp(ids);
    end
    MCSout=zeros(cfg.CCiter,2);
    
    for i = 1:cfg.CCiter
        MCSout(i,:)=estimate_max_clust(data,surro(:,i),cfg.type,inmask,[size(nii.img,1),size(nii.img,2),size(nii.img,3)],rCFT,i);
        save debug_this MCSout results
        disp([num2str(MCSout(1,1)) ' ' num2str(prctile(MCSout(2:i,1),95)) ' ' num2str(prctile(MCSout(2:i,2),95)) ' ' num2str(length(find(MCSout(2:i,1)>MCSout(1,1)))/length(MCSout(2:i,1)))])
        
    end
    
    results.MCSout=MCSout;
    results.rCFT=rCFT;
    results.CCT=prctile(MCSout(2:i,1),95);
    outclustcorr=zeros([size(nii.img,1),size(nii.img,2),size(nii.img,3)]);
    
    outclustcorr(inmask)=corr(data(:,inmask),surro(:,1),'type',cfg.type);
    
    tempmask=double(outclustcorr>=rCFT);

    a=bwlabeln(tempmask,18);
    tempU=unique(a);
    tempU(1)=[];
    rCCvol=zeros(size(a));
    for cs=1:length(tempU)
        temp=length(find(a==tempU(cs)));
        if(temp>=results.CCT)
            rCCvol(find(a==tempU(cs)))=1;
        end
    end
    
    results.rCCvol=rCCvol;
    results.a=a;
end



end

function out=estimate_max_clust(data,surro,type,inmask,kkk,rCFT,i)
disp(['Cluster correction iteration: ' num2str(i)])
surroR=zeros(kkk);
surroR(inmask)=corr(data(:,inmask),surro,'type',type);

tempmask=double(surroR>=rCFT);

a=bwlabeln(tempmask,18);
tempU=unique(a);
tempU(1)=[];
if(isempty(tempU))
    out=0;
else
    out=max(histc(a(:),tempU));
end
out(2)=max(surroR(inmask));


%     surro=zeros(iter,1);
%     parfor i=1:iter
% 		pe=randperm(size(xSim,1));
% 		temp=xSim(pe,pe);
% 		surro(i)=corr(temp(ids),modelMat(ids),'type',type);
%    	end
% 	[fi xi]=ksdensity(surro,'function','cdf','npoints',200);
% 	pval_left=interp1([-1 xi 1],[0 fi 1],out);    % trick to avoid NaNs
% 	pval=1-pval_left;
end



