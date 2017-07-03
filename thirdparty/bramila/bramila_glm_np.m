function cfg=bramila_glm_np(cfg)
% simple permutation based univariate parametric GLM for first level analysis
% Usage:
%   cfg = bramila_glm_np(cfg);
%
%   Input:
%   cfg.infile= path to nifti file
%   cfg.regressor = a 2D design matrix with number of rows equal to number of time
%       points, and number of columns equal to number of regressors. No
%       orthogonalization happens, i.e. regressors are tested independently
%   cfg.cdtP = p value of cluster defining threshold (default 0.05, uses t-distribution )
%   cft.cdtR = r value of cluster defining threshold (overrides cdtP)
%   cfg.NPERM = amount of permutations (default is 5000)
%   cfg.seed = random generator seed
%
%   Output:
%   cfg.vol = one volume for each regressor, cluster corrected
%   cfg.cth = cluster threshold

% input testing - still missing



T=size(cfg.regressor,1);
N=size(cfg.regressor,2);
nii=load_nii(cfg.infile);
mask=sign(var(nii.img,0,4));
inmask=find(mask>0);

% estimates the average autocorrelation of regressors
for r=1:N
    DF(r,1)=bramila_autocorr(cfg.regressor(:,r),cfg.regressor(:,r));
end

AC=T./DF;
mAC=mean(AC);

cfg.AC=AC;
cfg.mAC=mAC;

block_size=round(mAC);

% estimate cfg.cdtR
rho=0:0.001:0.5;
for rhoID=1:length(rho);
    tempp(rhoID)=pvalPearson('r', rho(rhoID), DF); % it's not the mean(AC) but the DF
end
cfg.cdtR=rho(min(find(tempp<cfg.cdtP)))



% creating permuted regressors
surrog=zeros(T,cfg.NPERM+1,N);


rng(cfg.seed);
for r=1:N
    for perms=1:cfg.NPERM+1
        if(perms==1)
            surrog(:,perms,r)=cfg.regressor(:,r);
        else
            surrog(:,perms,r)=bootres(cfg.regressor(:,r),0,block_size);
        end
    end
end

cfg.surrog=surrog;

csurr=zeros(cfg.NPERM+1,N);
data=reshape(nii.img,[],T);
data=data'; % Time in 1st dimension
for r=1:N
    for perms=1:cfg.NPERM+1
        disp(['Perm # ' num2str(perms) ' (' num2str(r) ')']);
        rsurrtemp=corr(data(:,inmask),squeeze(surrog(:,perms,r)));
        rsurr=zeros(size(mask));
        rsurr(inmask)=rsurrtemp;
        if(perms==1)
            vol=rsurr;
        end
        
        [x y z]=ind2sub(size(mask),find(rsurr>cfg.cdtR));
        a=[x y z];
        clust=spm_clusters(a',18);
        
        if(perms==1)
            volclust=clust;
            vola=a;
        end
        
        uu=unique(clust);
        uu(find(uu==0))=[];
        cs=histc(clust,uu);
        
        if(isempty(max(cs)))
            csurr(perms,r)=0;
        else
            csurr(perms,r)=max(cs);
        end
    end
    cfg.cth(r,1)=prctile(csurr(2:end,r),100*(1-cfg.cdtP));
    
    newmask=zeros(size(mask));
    
    for c = 1:length(volclust);
        cids=find(c==volclust);
        csize=length(cids);
        
        if(csize<cfg.cth(r))
            for cc=1:length(cids)
                newmask(vola(cids(cc),1),vola(cids(cc),2),vola(cids(cc),3))=0;
            end
        else
            for cc=1:length(cids)
                newmask(vola(cids(cc),1),vola(cids(cc),2),vola(cids(cc),3))=c;
            end
        end
    end
    
    cfg.vol(:,:,:,r)=vol;
    cfg.cmask(:,:,:,r)=newmask;
end

cfg.csurr=csurr;





%% helpers

function tsout=bootres(ts,type,MBS)
L=length(ts);
if(type==0)
    notok=1;
    tempids=(1:L)';
    while(notok==1)
        B=MBS+round(100*rand);	% block size
        R=ceil(L/B);
        temp=zeros(R*B,1);
        temp(1:L)=1:L;
        temp=circshift(temp,round(MBS/2)+round((length(temp)-round(MBS/2))*rand));
        mat=reshape(temp,B,R);
        mat=mat(:,randperm(R));
        temp=mat(:);
        
        tempids=temp(find(temp>0));
        if(min(abs(tempids-(1:L)'))>round(MBS/2))
            notok=0;
        end
    end
    tsout=ts(tempids);
else
    tsout=circshift(ts,round(MBS/2)+round((L-round(MBS/2))*rand));
end




