function SS = SumOfSqaures(M,opt)

if nargin==1
    opt=1;
end

if opt==1
    %M = M-repmat(mean(M),size(M,1),1);
    %SS = (sum(M.^2));
    mu = mean(M,1);
    dm = bsxfun(@minus,M, mu);
    SS = (sum(dm.^2));
else
    %M = M-repmat(nanmean(M),size(M,1),1);
    %SS = (nansum(M.^2));
    
    mu = nanmean(M,1);
    dm = bsxfun(@minus,M, mu);
    SS = (nansum(dm.^2));
end