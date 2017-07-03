function stats = bramila_two_sample_test(data,design,niter,nworker,use_t_vals,use_kernel_pvals,MAX_ARRAY_SIZE)
% BRAMILA_TWO_SAMPLE_TEST, a modified BRAMILA_TTEST2_NP for faster computations with the expence of
% increased memory consumption, vectorized permutations and skipped NaN checks.
% For small (non-integer) datasets BRAMILA_TTEST2_NP is recommended.
%
%   - Usage:
%   stats = bramila_two_sample_test(data,design,niter,nworker,use_t_vals,use_simple_pvals,MAX_ARRAY_SIZE)
%   - Input:
%   data = a matrix where each column is a subject and each row is a
%       data-point for example a voxel intensity in fMRI, a node level
%       value in a network, etc. NaN values will be ignored.
%   design = a row vector containing the numbers 1 and 2 for the two groups
%   niter = number of permutations (recommended 5000)
%   nworker = number of requested workers (if no existing pool found), default: local profile configuration (OPTIONAL)
%   use_t_vals = compute t-vals instead of simple mean difference, default: false (OPTIONAL)
%   use_kernel_pvals = compute kernel density p-values instead of simple ones, default: false (OPTIONAL)
%   MAX_ARRAY_SIZE = upper limit for the largest float array used by workers, if this limit is reached, 
%       work is divided in smaller chunks, default: 1e+9 (OPTIONAL)
%
%   - Output:
%   stats = a struct with the subfields
%       pvals = p-values for each datapoint; it returns in order the p-values
%       for the right tail and for the left tail
%       tvals = T-values OR mean differences for datapoint, positive tvals mean group 1 > group 2
%       use_t_vals = parameter setting [true/false]
%       use_simple_pvals = parameter setting [true/false]
%
% also check BRAMILA_TTEST2_NP for further details
%
% JanneK 2014-07-25 (janne.kauttonen@aalto.fi)
 
if nargin<7
    MAX_ARRAY_SIZE = 1e+9; % maximum array size that still fits worker memory
end

% create local cluster
try
    myCluster = gcp('nocreate');
    if isempty(myCluster)
        %delete(gcp)
        myCluster = parcluster('local');
        if nargin>3 && ~isempty(nworker)
            myCluster.NumWorkers=nworker;
        end
        parpool(myCluster);
    end
    N_workers = myCluster.NumWorkers;
catch err % old matlab?
    if ~matlabpool('size')
        if nargin>3 && ~isempty(nworker)
            eval(['matlabpool local ',num2str(nworker)]);
        else
            matlabpool local
        end
    end
    N_workers = matlabpool('size');
end

%% let's do some input validation first
if nnz(isnan(data))>0
    error('Your data has NaN''s, remove them and then try again')
end

if nargin<5 || (nargin>4 && isempty(use_t_vals))
   use_t_vals=false;
end
use_t_vals=logical(use_t_vals);

if nargin<6 || (nargin>5 && isempty(use_kernel_pvals))
   use_kernel_pvals=false;
end
use_kernel_pvals=logical(use_kernel_pvals);

if size(data,1)<1000
    use_kernel_pvals=true;
end

if use_t_vals
   if nnz(floor(data)~=data)==0
       warning('data contains only integers, using t-vals is not recommended!')
   end
end

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
if(niter<=0)
    disp('The variable niter should be a positive integer, function will continue assuming niter=5000')
    niter=5000;
end

nx = size(data(:,g1),2); 
ny = size(data(:,g2),2); % a scalar, => a scalar call to tinv

stats.tvals=tt_np(data,g1,g2,nx,ny,use_t_vals);

% we stop here if bad t-vals are presents
if nnz(isnan(stats.tvals))>0
    error('NaN t-vals found, check your data (repeating identical values?)')
end

NC = size(data,1); % number of comparisons
ND = size(data,2);

fprintf('Data size: %i samples, %i (%i+%i) features\n',NC,ND,nx,ny);
if NC>1e5 && use_kernel_pvals
   warning('Simple proportional p-val computational is recommended for large data!');
end

% create sets
M = NC*niter;
N_sets = min(NC,max(floor(M/MAX_ARRAY_SIZE),N_workers));
sets = mod(1:NC,N_sets)+1;
sets = sets(randperm(NC));

tic;
fprintf('Slicing data (%i sets)... ',N_sets);
for s=1:N_sets
    set{s}=(sets==s);
    dataset{s}=data(set{s},:);
    setsize(s)=nnz(set{s});
    pvals_set{s}=zeros(setsize(s),2);   
    tvals_set{s}=stats.tvals(set{s});
end
fprintf('done (%.1f sec)\n',toc); 

tic;
fprintf('Starting permutations\n');
parfor s=1:N_sets
    NC_set = setsize(s);
    outiter=zeros(NC_set,niter,'single');
    data_temp=dataset{s};
    
    for iter=1:niter
        perm=randperm(ND);
        data_temp=data_temp(:,perm);
        tval=tt_np(data_temp,g1,g2,nx,ny,use_t_vals);
        outiter(:,iter)=tval;
    end        
    
    for i=1:NC_set        
        pvals_set{s}(i,:) = tt_np_pval(tvals_set{s}(i),outiter(i,:),niter,use_kernel_pvals,use_t_vals);
    end    
    
    fprintf('...set nr. %i (of %i) done\n',s,N_sets);      
end
fprintf('done in %.1f sec\n',toc);    

tic;
fprintf('Pooling results... ');
stats.pvals = zeros(NC,2);
for s=1:N_sets
    stats.pvals(set{s},:)=pvals_set{s};
end
fprintf('done (%.1f sec)\n',toc); 

if nnz(isnan(stats.pvals))>0
    warning('NaN''s found in final p-vals')
end

stats.use_t_vals=use_t_vals;
stats.use_kernel_pvals=use_kernel_pvals;

end

function tval=tt_np(data,g1,g2,nx,ny,use_t_vals)
    % helper function similar to matlab function ttest2.m for the case of
    % groups with difference variance

    difference = mean(data(:,g1),2) - mean(data(:,g2),2);
    
    if use_t_vals
        s2x = var(data(:,g1),[],2);
        s2y = var(data(:,g2),[],2);
        s2xbar = s2x ./ nx;
        s2ybar = s2y ./ ny;
        se = sqrt(s2xbar + s2ybar);
        tval = difference ./ se;
    else
        tval = difference;
    end

end

function pval = tt_np_pval(tval,outiter,niter,use_kernel_pvals,use_t_vals)

if use_t_vals && use_kernel_pvals
    NCDF=200;
    if(niter>5000)
        NCDF=round(200*niter/5000);
    end
    [fi,xi]=ksdensity(outiter,'function','cdf','npoints',NCDF);

    edge_L = atanh(-1+eps);
    edge_R = atanh(1-eps);
    if tval<=edge_L
        pval=[1,0];
    elseif tval>=edge_R
        pval=[0,1];
    else
        pval_left=interp1([edge_L,xi,edge_R],[0,fi,1],tval);
        pval_right=1-pval_left;
        pval=[pval_right pval_left];
    end
elseif ~use_t_vals && use_kernel_pvals
   NCDF=200;
    if(niter>5000)
        NCDF=round(200*niter/5000);
    end
    [fi,xi]=ksdensity(outiter,'function','cdf','npoints',NCDF);
        
    edge_L = xi(1);
    edge_R = xi(end);
    if tval<=edge_L
        pval=[1,0];
    elseif tval>=edge_R
        pval=[0,1];
    else
        pval_left=interp1(xi,fi,tval);
        pval_right=1-pval_left;
        pval=[pval_right,pval_left];
    end
else    
    a=nnz(outiter>tval)/niter;
    pval=[a,1-a];    
end

end