function bspm_optimalhpf(r, filterRange)
% OPTIMALHPF
%
% This is a tool for assessing how different high-pass filters affect the
% estimation efficiency of an SPM design matrix. This function MUST be run
% in a folder containing an SPM.mat design file. In addition, this calls
% SPM functions. 
%
% USER SPECIFIED VARIABLES (below)
%   contrasts = contrasts of interest (each row represents one contrast)
%   conWeights = weights corresponding to importance of each contrast in 'contrasts'
%
% USAGE: optimalhpf(r, filterRange)
%   r = which run to estimate on? (default = 1)
%   filterRange = range of filters to assess (default = 5:500)
%
% EXAMPLE
%   In this example, there are three contrasts of interest (semi-colon
%   starts new row). The first contrast is weighted twice as much as the
%   second and third contrasts (because it is more important). It will be
%   executed for RUN 2 and for the filter range 100:200. 
%   
%       contrasts = [1 -1; 1 0; 0 1];
%       conWeights = [1 .5 .5];
%       >> optimalhpf(2, 100:200)
% 
%
% Bob Spunt, UCLA
% June 1, 2012 -- Created
% June 4, 2012 -- Added to FUNC + Presented in SCAN Lab Workshop

% ------------------------------- Copyright (C) 2014 -------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014


%---------USER INPUTS---------%
contrasts=[1 0 0 -1 0 0; 0 1 0 0 -1 0; 0 0 1 0 0 -1; -.5 .5 0 .5 -.5 0; -.5 0 .5 .5 0 -.5]';             % contrasts of interest (each row = new contrast)
conWeights=[1 1 1 1 1];          % how to weight each contrast in overall efficiency estimation? 
%--------END USER INPUT--------%

% Check argument
%------------------------------------------------------------------------
if nargin<2
    r=1;
    filterRange=5:500;
end
        
% Load design matrix
%------------------------------------------------------------------------
ncontrasts=size(contrasts,2);     
load SPM.mat
numruns = length(SPM.Sess);
if numruns==1 && r>1
    fprintf('\nWARNING! Found only 1 run in SPM.mat. Will run on run 1.\n');
end
K.RT = SPM.xY.RT;
Xmatrix = SPM.xX.X(SPM.Sess(r).row,SPM.Sess(r).col);
rundur = K.RT*length(Xmatrix);

% Begin filter loop
%------------------------------------------------------------------------
for i = 1:length(filterRange)

    
    K.HParam = filterRange(i);
    K.row = 1:length(Xmatrix);
    K = spm_filter(K);
    FXmatrix = spm_filter(K,Xmatrix);

    % compute efficiency for this iteration
    desmtx = FXmatrix;
    desmtx(:,end+1)=1;  % add intercept to desmtx
    ncols = size(desmtx,2);
    L = zeros(ncols,ncontrasts);
    L(1:size(contrasts,1),:) = contrasts;
    
    for c=1:ncontrasts
        
       efficiency(c)=1/trace(L(:,c)'*pinv(desmtx'*desmtx)*L(:,c));
       
    end
    
    effTotal=sum(conWeights.*efficiency);
    
    results(i,1) = filterRange(i);
    results(i,2) = effTotal;

end;
   
% Plot the results
%------------------------------------------------------------------------
figure('Color','white')
subplot(1,3,1); imagesc(zscore(Xmatrix)); colormap('gray');
title(sprintf('Run duration is %d s', rundur));
subplot(1,3,2); plot(results(:,1),results(:,2)); 
xlabel('High-Pass Filter (s)'); ylabel('Estimation Efficiency');
subplot(1,3,3); plot(results(2:end,1),diff(results(:,2))); 
xlabel('High-Pass Filter (s)'); ylabel('Gain in Estimation Efficiency');
large = find(diff(results(:,2))==max(diff(results(:,2)))) + 1;
title(sprintf('Largest gain is at HPF of %d s', results(large,1)));




 
 
 
 
