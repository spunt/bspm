function [clean_vol,r2] = bramila_regress(cfg)

% INPUT:
%   cfg.vol = 4D volume
%   cfg.reg = regressors time series (time in 1st dimension)
%   cfg.mask = 3D mask
% OUTPUT:
%   clean_vol = 4D volume after removal of regressors
%   r2 = correlation between original and fitted time-series

if ~isempty(cfg.reg)
    
    X_nuisance = cfg.reg;
    voldata = double(cfg.vol);
    
    X_nuisance=removezeros(X_nuisance);
    
    % for generality, add a contant term if not already present
    if ~hasconstant(X_nuisance)
        X_nuisance = zscore(X_nuisance);
        X_nuisance = [ones(size(X_nuisance,1),1),X_nuisance];
    end
    
    X_nuisance_pinv = pinv(X_nuisance);
   
    % prepare the 4D data to be regressed
    siz=size(voldata);
    tsdata=reshape(voldata,prod(siz(1:3)),siz(4))';

    T = size(tsdata,1);
    R = rank(X_nuisance);
    
    if T-R<=0
        error('!! Too few timepoints, all signal will be removed !! ( nuisance matrix rank >= number of timepoints )');
    elseif T/R<4
        warning('Number timepoints is small, lots of signal will be removed! ( timepoints/rank < 4 )');
    end
    
    m=mean(tsdata,1);	% storing the mean
    for row = 1:T
        tsdata(row,:)=tsdata(row,:)-m;
    end    

    if isfield(cfg,'mask') && all(size(cfg.mask)==siz(1:3))
        mask=cfg.mask;
    else
        mask = true(siz(1:3));
    end
   
    maskID=find(mask>0);
    N=length(maskID);
    if N~=nnz(mask) % sanity check
       error('Number of mask non-zeros mismatch!') ;
    end

    tsdataout=zeros(size(tsdata));
    r2 = zeros(1,size(tsdataout,2));

    fprintf('..Regression progress (%ik voxels, %i regressors): ',round(N/1000),size(X_nuisance,2)-1);    
    
    printiter=round(linspace(N/10,N*9/10,10));
    kk=1;
    for k=1:N    
        if kk<10 && k==printiter(kk)
            fprintf('%i%% ',kk*10);
            kk=kk+1;
        end

        y=tsdata(:,maskID(k));
        fit=X_nuisance*(X_nuisance_pinv*y);
        tsdataout(:,maskID(k))=y-fit;

        fit = fit - mean(fit);
        
        % compute correlation
        r = sum(sum(y.*fit))/sqrt(sum(sum(y.*y))*sum(sum(fit.*fit)));
        r2(maskID(k))=r^2;
        
    end

    %tsdataout=tsdataout+m;
    for row = 1:T
        tsdataout(row,:)=tsdataout(row,:)+m;
    end           
    tsdataout=tsdataout';
    clean_vol=reshape(tsdataout,size(voldata));
    r2 = reshape(r2,siz(1:3));
    r2(isnan(r2))=0;    
    
    fprintf('100%%\n....min(r2)=%3.2f, max(r2)=%3.2f, mean(r2)=%3.2f, std(r2)=%3.2f\n',min(r2(maskID)),max(r2(maskID)),mean(r2(maskID)),std(r2(maskID)));
    if max(r2(maskID))==1.0
       warning('A perfect match between nuisance regressors and %i voxel timeseries found, check your parameters !!',nnz(r2(maskID)==1));
    end

else
    
    fprintf('..No regressors, skipping regression\n')
    clean_vol = cfg.vol;
    s=size(clean_vol);
    r2=zeros(s(1),s(2),s(3));
    
end

end

function res = hasconstant(mat)

for i=1:size(mat,2)
    a=unique(mat(:,i));
    if length(a)==1 && a(1)~=0
        res = 1;
        return;
    end
end
res = 0;
end

function mat = removezeros(mat)

nullcols = [];
for i=1:size(mat,2)
    a=unique(mat(:,i));
    if length(a)==1 && a(1)==0
        nullcols(end+1)=i;
    end
end
mat(:,nullcols)=[];

if ~isempty(nullcols)
    fprintf('removed %i all-zero columns\n',length(nullcols));
end
end
