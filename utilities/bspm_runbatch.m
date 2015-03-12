function bspm_runbatch(job)
% BSPM_RUNBATCH
%
%   USAGE: bspm_runbatch(job)
%
spm_jobman('initcfg'); 
spm_jobman('run',job);
end

 
 
 
 
