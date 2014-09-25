function bspm_level1_parfor(jobs, ncores)
% BSPM_LEVEL1_PARFOR
%
%   ARGUMENTS:
%
%       jobs:   jobs structure
%       ncores: n cores to use
%

% ----------- Copyright (C) 2014 -----------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<2, error('bspm_level1_parfor(jobs, ncores)'); end

try
    spm('defaults','fmri'); spm_jobman('initcfg');
    matlabpool('open', ncores);
    parfor i = 1:length(jobs)
        fprintf('\nWorking on %s...', jobs(i).name);
        bspm_level1(jobs(i).images, jobs(i).general_info, jobs(i).runs, jobs(i).contrasts)
        fprintf('DONE\n'); 
    end
    matlabpool('close');
catch
    disp('SOMETHING BAD HAPPENED!');
    matlabpool('close');
    return
end
 
 
 
 
