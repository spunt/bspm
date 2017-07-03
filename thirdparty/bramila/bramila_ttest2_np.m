function stats = bramila_ttest2_np(data,design,niter)
% BRAMILA_TTEST2_NP - A "non parametric" two-sample T-test that, instead of
% relying on the t-distribtion, uses permutations of group labels to
% estimate the null distribution. The null distribution is computed
% independently for each data point (= row), i.e. we do not assume the same
% distribution for each datapoint. However, we do assume that the data
% points are comparable (e.g. they correspond to the same location
% collected across all subjects)
%
%   - Usage:
%   stats = bramila_ttest2_np(data,design,niter)
%   - Input:
%   data = a matrix where each column is a subject and each row is a
%       data-point for example a voxel intensity in fMRI, a node level
%       value in a network, etc. NaN values will be ignored.
%   design = a row vector containing the numbers 1 and 2 for the two groups
%   niter = number of permutations (recommended 5000)
%
%   - Output:
%   stats = a struct with the subfields
%       pvals = p-values for each datapoint; it returns in order the p-values
%       for the right tail and for the left tail
%       tvals = T-values for datapoint, positive tvals mean group 1 > group 2
%
% Notes: the null distribution is estimated using the matlab function
% ksdensity by interpolating the permuted data. The distribution is
% estimated over 200 points if niter<=5000, otherwise it is estimated over
% round(200*niter/5000) points, for greater precision.

% EG 2014-01-13 - enrico.glerean@aalto.fi

%% let's do some input validation first
Nsubj = size(data,2);   % number of subjects
if(size(design,2) ~= Nsubj)
    error('Mismatched number of subjects: the number of columns of data variable  should match the number of columns of the design variable')
end
if(size(design,1) ~= 1)
    error('The design variable should only contain 1 row')
end
% let's get the group IDs
g1 = find(design==1);
g2 = find(design==2);
if((length(g1)+length(g2))~=Nsubj)
    error('The design variable should only contain numbers 1 and 2')
end
% a check on niter
if(niter<0)
    disp('The variable niter should be a positive integer, function will continue assuming niter=5000')
    niter=5000;
end
if(niter==0)
    disp('I will not do permutations');
end
stats.tvals=tt_np(data,g1,g2);
% computing pvalues
NC = size(data,1); % number of comparisons
pvals=zeros(NC,2);
if(niter>0)
parfor n=1:NC % we treat each comparison independently
    pvals(n,:) = tt_np_pval(data(n,:),g1,g2,niter,stats.tvals(n));
end
end
stats.pvals=pvals;
end

function tval=tt_np(data,g1,g2)
    % helper function similar to matlab function ttest2.m for the case of
    % groups with difference variance
    xnans = isnan(data(:,g1));
    if any(xnans(:))
        nx = sum(~xnans,2);
    else
        nx = size(data(:,g1),2); 
    end
    ynans = isnan(data(:,g2));
    if any(ynans(:))
        ny = sum(~ynans,2);
    else
        ny = size(data(:,g2),2); % a scalar, => a scalar call to tinv
    end

    difference = nanmean(data(:,g1),2) - nanmean(data(:,g2),2);
    
    s2x = nanvar(data(:,g1),[],2);
    s2y = nanvar(data(:,g2),[],2);
    s2xbar = s2x ./ nx;
    s2ybar = s2y ./ ny;
    se = sqrt(s2xbar + s2ybar);
    if(any(se == 0) || any(isnan(se)))
        se(find(se==0))=Inf;
		se(find(isnan(se)))=Inf;
		%error('Group variance seems to be null or NaN, please check your data')
    end
    tval = difference ./ se;

end

function pval = tt_np_pval(data,g1,g2,niter,tval)
    outiter=zeros(niter,1);
    ND=length(data);
    for iter=1:niter
        perm=randperm(ND);
        % one could add a test to see that they are indeed permuted
        temp=data(perm);
        outiter(iter)=tt_np(temp,g1,g2);
    end
    NCDF=200;
    if(niter>5000)
        NCDF=round(200*niter/5000);
    end
    [fi xi]=ksdensity(outiter,'function','cdf','npoints',NCDF);
    
    % trick to avoid NaNs, we approximate the domain of the CDF between
    % -Inf and Inf using the atanh function and the eps matlab precision
    % variable
    
    pval_left=interp1([atanh(-1+eps) xi atanh(1-eps)],[0 fi 1],tval); 
    pval_right=1-pval_left;
    pval=[pval_right pval_left];
end


