function out = spm_rwls_run_plot(job)
% Does Plotting of movement parameter and residual statistic on normal GLMS 
% For outlier detection. 
% If model is estimated under wls, then the plot is extended to plot varince parameter and rewighted residuals after 
% the reweihting. 
% Model needs to be estimated under rwls toolbox 
% ______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_rwls_run_plot.m 3327 2010-05-18 08:27:32Z joern $



%-Load SPM.mat file
%-----------------------------------------------------------------------
SPM = [];
load(job.spmmat_plot{:});

spm_rwls_resstats(SPM,job.plot_subset,job.movparam); 
